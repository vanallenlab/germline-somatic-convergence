#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-Present Samantha Hoffman, Ryan Collins, and the Van Allen Lab @ DFCI  
# Distributed under terms of the GPL-2.0 License (see LICENSE)

"""
Identify pairs of related genes from two input gene lists
"""

import argparse
import itertools
import pandas as pd
from numpy import nanmin
from sys import stdout


# Declare constants
ppi_tiers = {'direct_functional' : 'direct_ppi_known_function', 
             'direct_unknown' : 'direct_ppi_unspecified_function', 
             'indirect' : 'nonspecific_ppi'}
criteria_tiers = {'same_gene' : 1,
                  'ligand_receptor' : 2,
                  'direct_ppi_known_function' : 2,
                  'direct_ppi_unspecified_function' : 3,
                  'nonspecific_ppi' : 4,
                  'protein_complex' : 4}


def load_genes(glist_in):
    """
    Load a list of genes
    """

    with open(glist_in) as fin:
        genes = set([g.rstrip() for g in fin.readlines()])

    return genes


def find_direct_overlaps(gg, sg):
    """
    Find pairs of genes that are present in both --germline and --somatic
    """

    return [(gene, gene, ) for gene in gg.intersection(sg)]


def find_ligand_receptor_pairs(gg, sg, lr_db):
    """
    Find pairs of genes that correspond to known ligand-receptor interactions
    """

    pairs = []

    # Find cases where germline gene is ligand and somatic gene is receptor
    gl_sr = (lr_db.ligand.isin(gg)) & (lr_db.receptor.isin(sg))
    pairs += (lr_db.loc[gl_sr, 'ligand receptor'.split()].values.tolist())

    # Find cases where somatic gene is ligand and germline gene is receptor
    sl_gr = (lr_db.ligand.isin(sg)) & (lr_db.receptor.isin(gg))
    pairs += (lr_db.loc[sl_gr, 'receptor ligand'.split()].values.tolist())

    return [(g, s, ) for g, s in pairs]


def find_complex_pairs(gg, sg, complexes):
    """
    Find germline & somatic pairs present in the same protein complex
    """

    pairs = []

    for cpx in complexes:
        g_hits = gg.intersection(cpx)
        s_hits = sg.intersection(cpx)
        for g in g_hits:
            for s in s_hits:
                if g != s:
                    pairs.append((g, s, ))

    return pairs


def format_output(identical_pairs, lr_pairs, ppi_pairs, complex_pairs, outfile,
                  report_counts=False):
    """
    Collapse a table of germline:somatic pairs and write as .tsv to outfile
    Optionally, if report_counts is True, returns a headerless row of counts 
    by criteria (order: any, same_gene, ligand_receptor, known_ppi, protein_complex)
    """

    # Dictionary for integrating results keyed on (germline, somatic)
    res = {}

    def _update_pair(res, pair, criteria):
        if pair not in res.keys():
            res[pair] = {'germline_gene' : pair[0],
                         'somatic_gene' : pair[1],
                         'criteria' : set([criteria])}
        else:
            res[pair]['criteria'].add(criteria)

        return res
    
    # Update results with identical genes found in both somatic & germline lists
    for pair in identical_pairs:
        res = _update_pair(res, pair, 'same_gene')

    # Update results with ligand-receptor interactions
    for pair in lr_pairs:
        res = _update_pair(res, pair, 'ligand_receptor')

    # Update results with known protein-protein interactions
    for ppi_tier, tier_name in ppi_tiers.items():
        for pair in ppi_pairs[ppi_tier]:
            res = _update_pair(res, pair, tier_name)

    # Update results with protein complex memberships
    for pair in complex_pairs:
        res = _update_pair(res, pair, 'protein_complex')

    # Annotate each pair with best tier
    for pair in res.keys():
        best_tier = nanmin([criteria_tiers[k] for k in res[pair]['criteria']])
        res[pair]['tier'] = best_tier

    # Format results as dataframe
    res_df = pd.DataFrame.from_dict(res, orient='index').reset_index(drop=True)
    if len(res_df) < 1:
        res_df = pd.DataFrame(columns='germline_gene somatic_gene criteria tier'.split())
    res_df.sort_values(by='tier germline_gene somatic_gene'.split(), axis=0, inplace=True)
    
    # Write dataframe to outfile
    if report_counts:
        counts_df = pd.DataFrame([(res_df.tier == 1).sum(), 
                                  (res_df.tier == 2).sum(), 
                                  (res_df.tier == 3).sum(),
                                  (res_df.tier == 4).sum()]).T
        counts_df.to_csv(outfile, sep='\t', index=False, header=False, na_rep='NA')
    else:
        if len(res_df) > 0:
            res_df.loc[:, 'criteria'] = res_df.criteria.apply(lambda x: ';'.join(sorted(list(x))))
            res_df.rename(columns = {res_df.columns[0] : '#' + res_df.columns[0]}, inplace=True)
        res_df.to_csv(outfile, sep='\t', index=False, na_rep='NA')


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--germline', help='list of germline genes', metavar='txt')
    parser.add_argument('--somatic', help='list of somatic genes', metavar='txt')
    parser.add_argument('--cellchat-db', help='curated .csv of Cellchat known ' +
                        'ligand-receptor pairs', metavar='txt')
    parser.add_argument('--ppi-db', help='curated .tsv of known protein-' +
                        'protein interactions as generated by curate_ebi_ppis.py', 
                        metavar='txt')
    parser.add_argument('--protein-complexes', help='List of known protein ' +
                        'complexes as generated by curate_wbi_complexes.py',
                         metavar='tsv')
    parser.add_argument('--out-tsv', help='output .tsv', metavar='tsv', 
                        default='stdout')
    parser.add_argument('--report-counts', default=False, action='store_true',
                        help='report counts of overlapping gene pairs per ' +
                        'criterion [default: return verbose list of all ' + 
                        'overlapping pairs]')
    args = parser.parse_args()

    # Load gene lists
    gg = load_genes(args.germline)
    sg = load_genes(args.somatic)

    # 1. Find genes present in both gg and sg
    identical_pairs = find_direct_overlaps(gg, sg)

    # 2. Find receptor-ligand pairs
    if args.cellchat_db is not None:
        lr_db = pd.read_csv(args.cellchat_db, sep=',')
        lr_pairs = find_ligand_receptor_pairs(gg, sg, lr_db)
    else:
        lr_pairs = []

    # 3. Find pairs with known protein-protein interactions
    # As of Dec 16, 2024, we now divide PPIs into three tiers of support
    if args.ppi_db is not None:
        ppi_df = pd.read_csv(args.ppi_db, sep='\t')
        ppi_pairs = {}
        for ppi_tier in ppi_tiers.keys():
            ppi_members_sub = ppi_df.loc[ppi_df['type'] == ppi_tier, 'members']
            ppis = [set(x.split(';')) for x in ppi_members_sub.values.tolist()]
            ppi_pairs[ppi_tier] = find_complex_pairs(gg, sg, ppis)
    else:
        ppi_pairs = {k : [] for k in ppi_tiers.keys()}

    # 4. Find pairs involved in the same protein complex
    if args.protein_complexes is not None:
        complex_df = pd.read_csv(args.protein_complexes, sep='\t')
        complexes = [set(x.split(';')) for x in complex_df.iloc[:, 1].values.tolist()]
        complex_pairs = find_complex_pairs(gg, sg, complexes)
    else:
        complex_pairs = []

    # Combine the results from 1-4 and output to --out-tsv
    if args.out_tsv in 'stdout /dev/stdout -'.split():
        outfile = stdout
    else:
        outfile = args.out_tsv
    format_output(identical_pairs, lr_pairs, ppi_pairs, complex_pairs, outfile, 
                  args.report_counts)


if __name__ == '__main__':
    main()

