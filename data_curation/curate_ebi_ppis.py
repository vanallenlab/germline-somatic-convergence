#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-Present Ryan Collins and the Van Allen Lab @ DFCI  
# Distributed under terms of the GPL-2.0 License (see LICENSE)

"""
Parse EBI IntAct database XML to extract lists of known protein-protein interactions
"""

import argparse
import pandas as pd
import xml.etree.ElementTree as ET
from os.path import basename
from re import sub

# Declare list of non-specific interaction types to exclude during curation
nonspec_itypes = ['association', 'colocalization', 'proximity', 
                  'self interaction', 'putative self interaction']


def get_ppis(xml_in, return_dataframe=False, eligible_genes=None):
    """
    Parse a single XML file and extract all PPIs asserted as "direct interaction"
    or "physical association"

    Returns a dict of {interaction_name : set(interactors)}
    Alternatively returns pd.DataFrame if return_dataframe is True
    """

    tree = ET.parse(xml_in)

    root = tree.getroot()

    res = {}

    # As discovered on December 5, 2024, each .xml can actually contain multiple
    # PPIs. These are stored in the <entry> element in .xml. We need to iterate 
    # over all of them and process each
    for i in range(len(root)):

        # First, build a dictionary mapping interactorRef codes to gene symbols
        # This is necessary because the actual PPIs held in the interactionList 
        # don't specify gene symbols, just references back to elements in interactorList
        interactors = {}
        ilist = [x for x in root[i].findall('./') if x.tag.endswith('interactorList')][0]

        # Get gene symbols for all interactors and store in interactors dict
        for interactor in ilist:

            inumber = interactor.attrib.get('id')

            symbols = set()

            inames = [e for e in interactor.findall('./') if e.tag.endswith('names')][0]
            for sub_e in inames.iter():
                if sub_e.tag.endswith('alias'):
                    if 'gene name' in sub_e.get('type'):
                        gene = sub_e.text.upper()
                        if ' ' not in gene:
                            symbols.add(gene)

            interactors[inumber] = symbols

        # The interactions themselves are stored in the interactionList element, 
        # which is a grandchild of the root
        ilist = [x for x in root[i].findall('./') if x.tag.endswith('interactionList')][0]

        # Each interaction is a single child "interaction" element of interactionList
        for interaction in ilist:

            # Check criteria for interaction and only keep direct interactions or
            # physical associations
            itype = [e for e in interaction.findall('./') if e.tag.endswith('interactionType')][0][0][0].text
            if itype in nonspec_itypes:
                continue

            # Get interaction name
            iname = [e for e in interaction.findall('./') if e.tag.endswith('names')][0][0].text

            # Get participants
            members = set()
            participant_list = [e for e in interaction.findall('./') if e.tag.endswith('participantList')][0]
            for participant in participant_list:
                # Get reference number to interactor
                try:
                    inumber = [x for x in participant.findall('./') if x.tag.endswith('interactorRef')][0].text
                    members.update(interactors[inumber])
                except:
                    # In at least one case, there was no interactorRef specified
                    # These situations seem to be pretty rare and should probably be skipped
                    continue

            # Filter members to eligible gene list, if optioned
            if eligible_genes is not None:
                members = members.intersection(eligible_genes)

            # Update entry in results aggregator
            if len(members) > 1:
                res[iname] = members

    # Convert res to pd.DataFrame if optioned
    if return_dataframe:
        for k, v in res.items():
            res[k] = ';'.join(sorted(list(v)))
        res = pd.DataFrame.from_dict(res, orient='index').reset_index()

    return res


def format_output_tsv(res):
    """
    Convert list of pd.DataFrame into two-column pd.DataFrame
    """

    # First, collapse all dataframes
    res_df = pd.concat(res, ignore_index=True)
    res_df.columns = '#interaction members'.split()

    # Deduplicate by members
    return res_df.loc[~res_df.members.duplicated(), :].reset_index(drop=True)


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

    # Process each .xml file
    res = []
    for xml_in in args.xml:
        cpx_id = sub('\.xml$', '', sub('_human', '', basename(xml_in)))
        if cpx_id.endswith('_negative'):
            continue
        res.append(get_ppis(xml_in, return_dataframe=True, eligible_genes=elig))

    # Exclude any interactions listed in the "negative" XML files
    # Per EBI documentation, these are interactions that have been refuted by
    # published evidence, and it sounds like they set a pretty high bar for this
    # TODO: need to implement this in the future

    # Format output as a two-column pd.DataFrame of interaction id & gene symbols
    # Before writing to --out-tsv
    res_df = format_output_tsv(res)
    res_df.to_csv(args.out_tsv, sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':
    main()

