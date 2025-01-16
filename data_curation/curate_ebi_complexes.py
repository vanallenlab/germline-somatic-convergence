#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-Present Samantha Hoffman, Ryan Collins, and the Van Allen Lab @ DFCI  
# Distributed under terms of the GPL-2.0 License (see LICENSE)

"""
Parse EBI Complex Portal XML to extract lists of genes involved in the same complex
"""

import argparse
import pandas as pd
import xmltodict
from os.path import basename
from re import sub


def get_complex_members(xml_in, eligible_genes=None):
    """
    Parse a single XML file and extract all human gene symbols
    """

    members = set()

    # Gene symbols are stored in the interactorList element, which is a 
    # grandchild of the root. Each member of the complex is a single child
    # element of interactorList
    with open(xml_in) as fin:
        ilist = xmltodict.parse(fin.read())['entrySet']['entry']['interactorList']['interactor']
        if isinstance(ilist, dict):
            ilist = [ilist]

    for member in ilist:

        symbols = set()
        
        alias_infos = member['names'].get('alias', [])
        if isinstance(alias_infos, dict):
            alias_infos = [alias_infos]
        for ainfo in alias_infos:
            gene = ainfo.get('#text')
            if ' ' in gene or any(c.islower() for c in gene):
                continue
            symbols.add(gene.upper())

        # If --eligible-genes is optioned, enforce that filter here
        if eligible_genes is not None:
            symbols = symbols.intersection(eligible_genes)

        if len(symbols) > 0:
            members.update(symbols)

    return members


def format_output_tsv(res):
    """
    Convert dict of complex : member mappings into two-column pd.DataFrame
    """

    # First, collapse all values from res dict as sorted semicolon-delimmed string
    complex_ids = [x for x in res.keys()]
    for cid in complex_ids:
        genes = res.get(cid, set())
        if len(genes) < 2:
            res.pop(cid)
        else:
            res[cid] = ';'.join(sorted(list(genes)))

    # Next, render as pd.DataFrame
    res_df = pd.DataFrame.from_dict(res, orient='index').reset_index()
    res_df.columns = '#complex members'.split()

    return res_df


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('xml', help='one or more input .xmls', metavar='xml', nargs='+')
    parser.add_argument('-g', '--eligible-genes', help='optional list of ' +
                        'eligible gene symbols to be included')
    parser.add_argument('--out-tsv', help='output .tsv', metavar='tsv')
    args = parser.parse_args()

    # Load list of eligible genes, if optioned
    if args.eligible_genes is not None:
        with open(args.eligible_genes) as fin:
            elig = set([g.rstrip() for g in fin.readlines()])
    else:
        elig = None

    # One row per complex
    res = {}
    for xml_in in args.xml:
        cpx_id = sub('\.xml$', '', sub('_human', '', basename(xml_in)))
        cpx_members = get_complex_members(xml_in, elig)
        if len(cpx_members) > 1:
            res[cpx_id] = cpx_members

    # Format output as a two-column pd.DataFrame of complex id & members
    # Before writing to --out-tsv
    res_df = format_output_tsv(res)
    res_df.to_csv(args.out_tsv, sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':
    main()

