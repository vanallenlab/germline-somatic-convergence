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


def preprocess_gtf(gtf_in, elig_genes=None):
    """
    Reads and preprocesses a GTF
    Returns filtered GTF as a pbt.BedTool
    """

    # Load GTF as pbt.BedTool and clean chromosome prefixes
    def _clean_chr(feat):
        feat.chrom = sub('^chr', '', str(feat.chrom))
        return feat
    gtf_bt = pbt.BedTool(gtf_in).each(_clean_chr)

    # Filter GTF to include only eligible genes, if optioned
    def _gene_filter(feat, elig_genes=None):
        gname = feat.attrs.get('gene_name')
        if gname is None:
            return False
        return gname in elig_genes
    gtf_bt = gtf_bt.filter(_gene_filter, elig_genes=elig_genes)

    # Return a copy (not generator) of gtf as pbt.BedTool
    return gtf_bt.saveas()


def precompute_hits(gwas_bt, gtf_bt):
    """
    Identify all hits between GWAS peaks and gene features
    Returns: dicts mapping GWAS catalog idx to most "severe"/high-impact overlap
             and the gene(s) affected by that feature overlap (or nearest gene
             body for intergenic variants)
    """
    
    # Determine hits between gwas_bt and gtf_bt and save as pd.DataFrame
    hits_df = gwas_bt.intersect(gtf_bt, wao=True).\
                      to_dataframe(disable_auto_names=True, header=None).\
                      iloc[:, [3, 6, 12]].\
                      rename(columns={3 : 'idx', 6 : 'feature', 12 : 'attribs'})

    # Extract gene symbol for every hit
    def _get_symbol(attribs):
        for ae in str(attribs).split(';'):
            if ae.startswith('gene_name'):
                return sub('"', '', ae.split(' ')[1])
        return None
    hits_df['gene'] = hits_df['attribs'].astype(str).apply(_get_symbol)
    hits_df.drop('attribs', axis=1, inplace=True)

    # Precompute dict mapping GWAS catalog index to nearest gene body
    gene_bt = gtf_bt.filter(lambda x: x.fields[2] in 'gene transcript'.split()).sort()
    prox_df = gwas_bt.sort().closest(gene_bt, t='all').\
                      to_dataframe(disable_auto_names=True, header=None).\
                      iloc[:, [3, 12]].\
                      rename(columns={3 : 'idx', 12 : 'attribs'})
    prox_df['gene'] = prox_df['attribs'].astype(str).apply(_get_symbol)

    # Build dicts keyed by GWAS catalog index mapping that index to the most
    # severe consequence and corresponding gene symbol(s) based on GTF feature overlap
    gene_context = {}
    mapped_genes = {}
    for idx in hits_df.idx.unique():
        idx_hits = hits_df.loc[hits_df.idx == idx, 'feature gene'.split()]
        worst = None
        for ftype in 'CDS coding_exon UTR other_exon'.split():
            if any(idx_hits.feature == ftype):
                worst = ftype
                genes = idx_hits.loc[idx_hits.feature == ftype, 'gene']
                gstr = ','.join(sorted(list(set(genes.values))))
                break
        if worst is None:
            if any(idx_hits.feature.isin('gene transcript'.split())):
                worst = 'intron'
                genes = idx_hits.loc[idx_hits.feature.isin('gene transcript'.split()), 'gene']
                gstr = ','.join(sorted(list(set(genes.values))))
            else:
                worst = 'intergenic'
                genes = prox_df.loc[prox_df.idx == idx, 'gene']
                gstr = ','.join(sorted(list(set(genes.values))))
        gene_context[idx] = worst
        mapped_genes[idx] = gstr

    return gene_context, mapped_genes


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
    parser.add_argument('--eligible-genes', metavar='txt', 
                        help='List of gene symbols to consider in --gtf')
    args = parser.parse_args()

    # Load input tsv and format as pbt.BedTool
    gwas_df, gwas_bt = load_raw_tsv(args.tsv_in)

    # Load eligible genes, if optioned
    if args.eligible_genes is not None:
        with open(args.eligible_genes) as fin:
            elig_genes = set([g.rstrip() for g in fin.readlines()])
    else:
        elig_genes = None

    # Load and preprocess GTF
    gtf_bt = preprocess_gtf(args.gtf, elig_genes)

    # Precompute overlap between GWAS catalog coordinates and --gtf
    gene_context, mapped_genes = precompute_hits(gwas_bt, gtf_bt)

    # Update input tsv with coding|noncoding distinction
    gwas_df.loc[:, 'gene_context'] = gwas_df.idx.map(gene_context)

    # Update MAPPED_GENE with annotations from --gtf
    gwas_df.loc[:, 'MAPPED_GENE'] = gwas_df.idx.map(mapped_genes)

    # Write annotated data to --tsv-out
    gwas_df.to_csv(args.tsv_out, index=False, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()

