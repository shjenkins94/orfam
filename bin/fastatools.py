#!/usr/bin/env python

import argparse

from Bio import SeqIO


def _translate(args):
    with open(args.fasta, "r") as in_fasta:
        for record in SeqIO.parse(in_fasta, 'fasta'):
            print(">" + record.id)
            print(record.seq.translate())


def _extract(args):
    with open(args.ids, "r") as in_ids:
        extract_ids = in_ids.read().splitlines()
    with open(args.fasta, "r") as in_fasta:
        for record in SeqIO.parse(in_fasta, 'fasta'):
            if record.id in extract_ids:
                print(record.format('fasta'))


def _exclude(args):
    with open(args.ids, "r") as in_ids:
        exclude_ids = in_ids.read().splitlines()
    with open(args.fasta, "r") as in_fasta:
        for record in SeqIO.parse(in_fasta, 'fasta'):
            if record.id not in exclude_ids:
                print(record.format('fasta'))


def _id_list(args):
    with open(args.fasta, "r") as in_fasta:
        for record in SeqIO.parse(in_fasta, 'fasta'):
            print(record.id)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="FASTA sequence manipulation tools")
    subparser = parser.add_subparsers(help="sub-commands")

    # parse translation arguments
    parser_translate = subparser.add_parser(
        'translate', help=("translate FASTA file"))
    parser_translate.add_argument(
        'fasta', help=("FASTA file"))
    parser_translate.set_defaults(func=_translate)

    parser_extract = subparser.add_parser(
        'extract', help=("extract sequence by ids"))
    parser_extract.add_argument(
        'fasta', help=("FASTA file"))
    parser_extract.add_argument(
        'ids', help=("one id in each line"))
    parser_extract.set_defaults(func=_extract)

    parser_exclude = subparser.add_parser(
        'exclude', help=("exclude sequence by ids"))
    parser_exclude.add_argument(
        'fasta', help=("FASTA file"))
    parser_exclude.add_argument(
        'ids', help=("ids to remove (txt)"))
    parser_exclude.set_defaults(func=_exclude)

    parser_id_list = subparser.add_parser(
        'id_list', help=("list sequence ids"))
    parser_id_list.add_argument(
        'fasta', help=("FASTA file"))
    parser_id_list.set_defaults(func=_id_list)

    # Parse arguments and run the sub-command
    args = parser.parse_args()
    args.func(args)
