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
import xmltodict
from os.path import basename
from re import sub

# Declare list of less specific interaction types
direct_nofunc_itypes = ['direct interaction', 'putative self interaction', 'self interaction']
indirect_itypes = ['physical association']
drop_itypes = ['colocalization', 'proximity', 'association']


def get_ppis(xml_in, return_dataframe=False, eligible_genes=None, max_members=10e10):
    """
    Parse a single XML file and curate all PPIs

    Returns a dict of {interaction_tier : {interaction_name : set(interactors)} }
    Alternatively returns pd.DataFrame if return_dataframe is True
    """

    # Instantiate empty dict for holding curated results
    res = {'direct_functional' : dict(), 'direct_unknown' : dict(), 
           'indirect' : dict()}

    # Read contents from XML into dict
    with open(xml_in) as fin:
        entries = xmltodict.parse(fin.read())['entrySet']['entry']
        if isinstance(entries, dict):
            entries = [entries]

    # As discovered on December 5, 2024, each .xml can actually contain multiple
    # PPIs. These are stored in the <entry> element in .xml. We need to iterate 
    # over all of them and process each
    for entry in entries:
        
        # First, build a dictionary mapping interactorRef codes to gene symbols
        # This is necessary because the actual PPIs held in the interactionList 
        # don't specify gene symbols, just references back to elements in interactorList

        interactors = dict()
        all_members = entry['interactorList']['interactor']
        if isinstance(all_members, dict):
            all_members = [all_members]
        
        # Always filter to human proteins only, and only continue if there are 
        # at least two interactors specified
        all_members = [m for m in all_members if \
                       m['organism']['names']['shortLabel'] == 'human' and \
                       m['interactorType']['names']['shortLabel'] == 'protein']
        if len(all_members) < 2:
            continue

        # Get gene symbols for all interactors and store in interactors dict
        for member in all_members:
            inumber = int(member['@id'])
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
                interactors[inumber] = symbols

        # The interactions themselves are stored in the interactionList element, 
        # which is a grandchild of the root. Each interaction is a single child 
        # "interaction" element of interactionList
        ilist = entry['interactionList']['interaction']
        if isinstance(ilist, dict):
            ilist = [ilist]
        for interaction in ilist:

            # Get interaction name
            iname = interaction['names']['shortLabel']

            # Assign interaction specificity tier
            # Never keep cellular colocalization or proximity, as they are too coarse
            itype = interaction['interactionType']['names']['shortLabel']
            if itype in drop_itypes:
                continue
            elif itype in indirect_itypes:
                itier = 'indirect'
            elif itype in direct_nofunc_itypes:
                itier = 'direct_unknown'
            else:
                itier = 'direct_functional'

            # Get participants, and only retain interaction if >1 participant
            imembers = set()
            all_participants = interaction['participantList']['participant']
            if isinstance(all_participants, dict):
                all_participants = [all_participants]
            if len(all_participants) < 2:
                continue

            for participant in all_participants:
                try:
                    inumber = int(participant['interactorRef'])
                except:
                    # In at least one case, there was no interactorRef specified and
                    # it was not clear why, nor could I find documentation online.
                    # These situations seem to be pretty rare and should be skipped
                    continue
                if inumber in interactors.keys():
                    imembers.update(interactors[inumber])

            # Filter members to eligible gene list, if optioned
            if eligible_genes is not None:
                imembers = imembers.intersection(eligible_genes)

            # Update entry in results aggregator
            if len(imembers) > 1 and len(imembers) <= max_members:
                res[itier][iname] = imembers

    # Convert res to pd.DataFrame if optioned
    if return_dataframe:
        subres_dfs = []
        for itype, subres in res.items():
            if len(subres) == 0:
                continue
            for k, v in subres.items():
                subres[k] = ';'.join(sorted(list(v)))
            subres_df = pd.DataFrame.from_dict(subres, orient='index').reset_index()
            subres_df.columns = 'interaction members'.split()
            subres_df['type'] = itype
            subres_dfs.append(subres_df['interaction type members'.split()])
        if len(subres_dfs) > 0:
            res = pd.concat(subres_dfs, ignore_index=True)
        else:
            res = pd.DataFrame(columns='interaction type members'.split())

    return res


def format_output_tsv(res):
    """
    Convert list of pd.DataFrame into two-column pd.DataFrame
    """

    # First, collapse all dataframes
    res_df = pd.concat(res, ignore_index=True)
    res_df.columns = '#interaction type members'.split()

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
    parser.add_argument('-m', '--max-members', default=10e10, type=int,
                        help='Maximum number of members to allow in a PPI')
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
        res.append(get_ppis(xml_in, return_dataframe=True, eligible_genes=elig, 
                            max_members=args.max_members))

    # Exclude any interactions listed in the "negative" XML files
    # Per EBI documentation, these are interactions that have been refuted by
    # published evidence, and it sounds like they set a pretty high bar for this
    # TODO: could implement this in the future

    # Format output as a three-column pd.DataFrame of interaction id, tier, 
    # and member gene symbols before writing to --out-tsv
    res_df = format_output_tsv(res)
    res_df.to_csv(args.out_tsv, sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':
    main()

