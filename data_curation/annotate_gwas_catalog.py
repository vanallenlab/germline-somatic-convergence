#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-Present Ryan Collins and the Van Allen Lab @ DFCI  
# Distributed under terms of the GPL-2.0 License (see LICENSE)

"""
Annotate GWAS catalog .tsv versus GTF
"""


import argparse
import pandas as pd
import pybedtools as pbt
from re import sub


def load_raw_tsv(tsv_in):
    """
    Load a .tsv of results from the GWAS catalog and format as a pbt.BedTool
    Returns: both the original data frame and the corresponding pbt.BedTool
    """

    # Read data
    df = pd.read_csv(tsv_in, sep='\t')
    df.insert(0, 'idx', df.index.values)
    df['locations'] = \
        df.apply(lambda x: ';'.join(['{}:{}'.format(chrom, int(pos)) \
                                     for chrom, pos \
                                     in zip(str(x.CHR_ID).split(';'), \
                                            str(x.CHR_POS).split(';'))]), axis=1)

    # Extract coordinates and reformat as a pbt.BedTool with GWAS catalog index as feature name
    def _format_record(data):
        fname = str(data['idx'])
        cstrs = []
        for csub in data.locations.split(';'):
            coords = csub.split(':')
            chrom = str(coords[0])
            start = str(coords[1])
            stop = str(int(coords[1]) + 1)
            cstrs.append('\t'.join([chrom, start, stop, fname]))
        return '\n'.join(cstrs)
    coord_strings = df.apply(_format_record, axis=1).tolist()
    bt = pbt.BedTool('\n'.join(coord_strings), from_string=True)

    # Return both the original dataframe (to be written out after annotation) and
    # the pbt.BedTool of coordinates (for intersecting vs. gene annotations)
    return df, bt


def precompute_hits(gwas_bt, gtf_in):
    """
    Identify all hits between GWAS peaks and gene features
    Returns: dict mapping GWAS catalog idx to most "severe"/high-impact overlap
    """

    # Load GTF as pbt.BedTool and clean chromosome prefixes
    def _clean_chr(feat):
        feat.chrom = sub('^chr', '', str(feat.chrom))
        return feat
    gtf_bt = pbt.BedTool(gtf_in).each(_clean_chr)
    
    # Determine hits between gwas_bt and gtf_bt and save as pd.DataFrame
    hits_df = gwas_bt.intersect(gtf_bt, wao=True).\
                      to_dataframe(disable_auto_names=True, header=None).\
                      iloc[:, [3, 6]].\
                      rename(columns={3 : 'idx', 6 : 'feature'})

    # Build a dict keyed by GWAS catalog index mapping that index to the most
    # severe consequence based on simple GTF feature overlap
    hits_dict = {}
    for idx in hits_df.idx.unique():
        idx_hits = set(hits_df.loc[hits_df.idx == idx, 'feature'].tolist())
        if 'CDS' in idx_hits:
            worst = 'coding_exon'
        elif 'UTR' in idx_hits:
            worst = 'UTR'
        elif 'exon' in idx_hits:
            worst = 'other_exon'
        elif 'gene' in idx_hits or 'transcript' in idx_hits:
            worst = 'intron'
        else:
            worst = 'intergenic'
        hits_dict[idx] = worst

    return hits_dict


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tsv-in', help='input .tsv', metavar='tsv')
    parser.add_argument('--tsv-out', help='output .tsv', metavar='tsv')
    parser.add_argument('--gtf', help='GTF for annotating CNVs', required=True)
    args = parser.parse_args()

    # Load input tsv and format as pbt.BedTool
    gwas_df, gwas_bt = load_raw_tsv(args.tsv_in)

    # Precompute overlap between GWAS catalog coordinates and --gtf
    hits = precompute_hits(gwas_bt, args.gtf)

    # Update input tsv with coding|noncoding distinction
    gwas_df.loc[:, 'gene_context'] = gwas_df.idx.map(hits)

    # Write annotated data to --tsv-out
    gwas_df.to_csv(args.tsv_out, index=False, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()

