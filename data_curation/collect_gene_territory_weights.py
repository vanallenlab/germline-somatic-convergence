#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-Present Samantha Hoffman, Ryan Collins, and the Van Allen Lab @ DFCI  
# Distributed under terms of the GPL-2.0 License (see LICENSE)

"""
Build map of gene territories for permutation weighting
"""

import argparse
import pandas as pd
import pybedtools as pbt
from numpy import nansum


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bed-in', help='input .csv', metavar='bed', 
                        required=True)
    parser.add_argument('--genes-list', help='list of eligible genes', 
                        metavar='txt', required=True)
    parser.add_argument('--tsv-out', help='output .tsv', metavar='tsv',
                        required=True)
    args = parser.parse_args()

    # Load BED of gene territories as pd.DataFrame
    data = pd.read_csv(args.bed_in, sep='\t', names='chrom start end gene'.split())

    # Count number of unique basepairs per gene
    res = {}
    with open(args.genes_list) as glin:
        for gstr in glin.readlines():
            gene = gstr.rstrip()
            # Initialize gene with 1 nucleotide to ensure all genes get non-zero weights
            if gene not in res.keys():
                res[gene] = 1
            for reg in pbt.BedTool.from_dataframe(data[data.gene == gene]).merge():
                res[gene] += len(reg)
    
    # Format output as tsv
    res_df = pd.DataFrame.from_dict(res, orient='index').reset_index()
    res_df.columns = '#gene bp'.split()
    res_df.sort_values('#gene').to_csv(args.tsv_out, sep='\t', index=False)


if __name__ == '__main__':
    main()

