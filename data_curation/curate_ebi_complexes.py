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
import xml.etree.ElementTree as ET
from os.path import basename
from re import sub


def get_complex_members(xml_in):
    """
    Parse a single XML file and extract all human gene symbols
    """

    tree = ET.parse(xml_in)

    root = tree.getroot()

    members = set()

    # Gene symbols are stored in the interactorList element, which is a 
    # grandchild of the root
    ilist = [x for x in root[0].findall('./') if x.tag.endswith('interactorList')][0]

    # Each complex interactor is a single child element of interactorList
    for interactor in ilist:

        # We want to extract the "names" element within each interactor
        names = [x for x in interactor.findall('./') if x.tag.endswith('names')][0]

        # Each names element can hvae multiple alias sub-elements
        # We want to extract text from each
        for sub_e in names.iter():
            if sub_e.tag.endswith('alias'):
                if 'gene name' in sub_e.get('type'):
                    gene = sub_e.text
                    if ' ' not in gene and gene.isupper():
                        members.add(gene)

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
    parser.add_argument('--out-tsv', help='output .tsv', metavar='tsv')
    args = parser.parse_args()

    # One row per complex
    res = {}
    for xml_in in args.xml:
        cpx_id = sub('\.xml$', '', sub('_human', '', basename(xml_in)))
        res[cpx_id] = get_complex_members(xml_in)

    # Format output as a two-column pd.DataFrame of complex id & members
    # Before writing to --out-tsv
    res_df = format_output_tsv(res)
    res_df.to_csv(args.out_tsv, sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':
    main()

