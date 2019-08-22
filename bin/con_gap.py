#!/usr/bin/env python

import argparse
from Bio import AlignIO


# 3. The list of conserved regions is sorted by start position. Then j, which
#    represents the the position in the original reference sequence is set to
#    -1 and reg, representing the next region to locate, is set to zero.
#    reg_max is set to the number of conserved regions, so that the function
#    can stop as soon as the last conserved region is located. A for loop then
#    goes through aln_seq with two iterators, i and x. i represents the
#    position in the alignment, and x represents the amino acid at position i.
#    j only increases if x is not a gap "-". When j is equal to the start
#    position of the current region, i is added on to con_reg as the
#    start position of the conserverd region in the alignment. When j is equal
#    to the stop position of the current region - 1 (the position of the last
#    amino acid in the region), i + 1 (the index that python and bed do not
#    include) is set as the stop position, and reg is incremented by 1. When
#    reg = reg_max ,the start and stop positions of all conserved regions have
#    been found, and the updated con_reg array is returned.
def find_regs(aln_seq, con_reg):
    con_reg.sort(key=lambda x: x[1])
    j = -1
    reg = 0
    reg_max = len(con_reg)
    for i, x in enumerate(aln_seq):
        if x != "-":
            j += 1
            if j == con_reg[reg][1]:
                con_reg[reg][3] = i
            elif j == con_reg[reg][2] - 1:
                con_reg[reg][4] = i + 1
                reg += 1
                if reg == reg_max:
                    break
    return con_reg


# 2. Takes in arguments and makes:
#    or_aln an AlignIO alignment
#    rem_id a list of strings
#    ref_bed an array of (Ref OR ID, region start, region end, -1, -1)
#    Then it checks that all regions have the same chrom value. If not, error
#    The function then attempts to get the index of the reference or. If it
#    can't, error.
#    The alignment sequence of the reference OR, and the list of conserved
#    regions are sent to find_regs.
def con_gap(args):
    with open(args.or_aln, "r") as in_or_aln:
        or_aln = AlignIO.read(in_or_aln, "fasta")
    con_reg = []
    ref_id = []

    with open(args.ref_bed, "r") as in_ref_bed:
        for line in in_ref_bed:
            cols = line.rstrip().split()
            assert(len(cols) == 3)
            # Ref_OR ID, Ref TM Start, Ref TM Stop
            con_reg.append([cols[0], int(cols[1]), int(cols[2]), -1, -1])
            ref_id.append(cols[0])

    if len(set(ref_id)) != 1:
        raise ValueError("Features in -b BED file must have same chrom value.")

    ref_id = ref_id[0]
    ref_or_index = -1
    for i in range(0, len(or_aln)):
        if or_aln[i].id == ref_id:
            ref_or_index = i
            break

    if ref_or_index == -1:
        raise NameError("-R sequence header must start with -b chrom value" +
                        " (" + ref_id + ")")

    con_reg = find_regs(or_aln[ref_or_index].seq, con_reg)
# 6. For each region, or_aln is trimmed so that it only contains the current
#    region. The sequence of each record in the trimmed alignment is gone
#    through. If the record has a gap at a position and the reference OR
#    sequence does not, then reg_gaps increases. If a sequence has more than
#    arg.gaps number of gaps in any conserved region, its id is printed.
#    If the cumulative flag is true, then a dictionary is created with entries
#    for each record, and the value of each entry increases for each gap found
#    after all regions have been evaluated, records with more than arg.gaps
#    gaps are printed.
    if args.cumulative:
        c_gap = {}
        for record in or_aln:
            c_gap[record.id] = 0
        print(c_gap)
    for reg in con_reg:
        reg_start = reg[3]
        reg_end = reg[4]
        reg_align = or_aln[:, reg_start:reg_end]
        reg_align_ref_or = reg_align[ref_or_index]
        for record in reg_align:
            reg_gaps = 0
            for i in range(0, len(record.seq)):
                if reg_align_ref_or.seq[i] != "-" and record.seq[i] == "-":
                    reg_gaps += 1
                    if args.cumulative:
                        c_gap[record.id] += 1
            if reg_gaps > args.gaps:
                print(record.id)
    if args.cumulative:
        for record in or_aln:
            if c_gap[record.id] > args.gaps:
                print(record.id)


# 1. takes in arguments (or_aln, ref_bed, rem_id, gaps) and passes
#    them to exclue_con_gap
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="exclude OR candidates with long gaps in conserved" +
        " regions")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument(
        '-A', '--or_aln',
        required=True,
        help=("Result of aligning potential ORs with references (FASTA)"))
    required.add_argument(
        '-B', '--ref_bed',
        required=True,
        help=("Conserved regions of one of the reference sequences (BED)"))
    required.add_argument(
        '-G', '--gaps',
        required=True,
        type=int,
        help=("Maximum nummber of gaps an intact OR can have (int)"))
    optional.add_argument(
        '-c', '--cumulative',
        action='store_true',
        help=("if the gap cutoff is for each region or cumulative"))

    args = parser.parse_args()
    con_gap(args)
