#!/usr/bin/env python

import argparse
from Bio import AlignIO


# 3. The list of TMs is sorted by start position. Then j, which represents the
#    the position in the original reference sequence is set to -1 and TM,
#    representing the number of TMs in the bed file, is set to zero.
#    A for loop then goes through aln_seq with two iterators, i and x. i
#    represents the position in the alignment, and x represents the amino acid
#    at position i.
#    j only increases if x is not a gap "-". When j is equal to the start
#    position of the current TM, i is recorded as the start position of the
#    TM in the alignment. When j is equal to the stop position of the current
#    TM - 1 (the position of the last amino acid in the TM), i + 1 (the index
#    that python and bed do not include) is set as the stop position, and TM
#    is incremented by 1. When TM = 7 ,the start and stop positions of all
#    7 TMs have been found, and the updated or_tms array is returned.
def find_tms(aln_seq, or_tms):
    or_tms.sort(key=lambda x: x[1])
    j = -1
    TM = 0
    for i, x in enumerate(aln_seq):
        if x != "-":
            j += 1
            if j == or_tms[TM][1]:
                or_tms[TM][3] = i
            elif j == or_tms[TM][2] - 1:
                or_tms[TM][4] = i + 1
                TM += 1
                if TM == 7:
                    break
    return or_tms


# 2. Takes in arguments and makes:
#    ref_or a SeqIO Seqrecord,
#    orf_aln an AlignIO alignment
#    in_ref_or_tm an array of (Ref OR ID, TM start, TM end, -1, -1)
#    rem_id an array
#    Then it checks that all TMs have the same chrom value. If not, error
#    The function then attempts to get the index of the reference or. If it
#    can't, error.
#    The alignment sequence of REF OR, and the list of TM regions are sent to
#    find_tms.
def exclude_tm_long_gap(args):
    with open(args.orf_aln, "r") as in_orf_aln:
        orf_aln = AlignIO.read(in_orf_aln, "fasta")
    with open(args.rem_id, "r") as in_rem_id:
        remove_names = in_rem_id.read().splitlines()
    or_tms = []
    ref_id = []
    with open(args.ref_or_tm, "r") as in_ref_or_tm:
        for line in in_ref_or_tm:
            cols = line.rstrip().split()
            assert(len(cols) == 3)
            # Ref_OR ID, Ref TM Start, Ref TM Stop
            or_tms.append([cols[0], int(cols[1]), int(cols[2]), -1, -1])
            ref_id.append(cols[0])

    if len(set(ref_id)) != 1:
        raise ValueError("Features in -b BED file must have same chrom value.")

    ref_id = ref_id[0]
    ref_or_index = -1
    for i in range(0, len(orf_aln)):
        if orf_aln[i].id == ref_id:
            ref_or_index = i
            break

    if ref_or_index == -1:
        raise NameError("-R sequence header must start with -b chrom value" +
                        " (" + ref_id + ")")

    or_tms = find_tms(orf_aln[ref_or_index].seq, or_tms)

# 6. For each TM, orf_aln is trimmed so that it only contains the current TM.
#    The sequence of each record in the trimmed alignment is gone through.
#    If the record has a gap at a position and the reference or sequence does
#    not, then tm_gaps increases. If a sequence has more than 4 gaps in any tm,
#    its id is added to remove names. After all TMs have been checked,
#    sequences not on the remove_names list are formatted as GFF entries and
#    returned
    for tm in or_tms:
        tm_start = tm[3]
        tm_end = tm[4]
        tm_align = orf_aln[:, tm_start:tm_end]
        tm_align_ref_or = tm_align[ref_or_index]
        for record in tm_align:
            tm_gaps = 0
            for i in range(0, len(record.seq)):
                if tm_align_ref_or.seq[i] != "-" and record.seq[i] == "-":
                    tm_gaps += 1
            if tm_gaps > 4:
                remove_names.append(record.id)

    for record in orf_aln:
        if record.id not in remove_names:
            chrom = record.id.split(":")[0]
            ref_start = int(
                record.id.split(":")[1].split("(")[0].split("-")[0]) + 1
            ref_end = int(record.id.split(":")[1].split("(")[0].split("-")[1])
            strand = record.id.split("(")[1].split(")")[0]
            # output (GFF format)
            orf = [chrom, "tm_gap", "ORF",
                   str(ref_start), str(ref_end), "0.0", strand, ".", "."]
            print('\t'.join(orf))


# 1. takes in arguments (ref_or, orf_aln, ref_or_tm, and rem_id) and passes
#    them to exclue_tm_long_gap
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="exclude OR candidates with long gaps in transmembrane" +
        "regions")
    parser.add_argument(
        '--orf_aln',
        required=True,
        help=("Open Reading Frame (Results of multiple sequence alignment)"))
    parser.add_argument(
        '--ref_or_tm',
        required=True,
        help=("TMs of OR reference sequence (BED format)"))
    parser.add_argument(
        '--rem_id',
        required=True,
        help=("ids of reference or sequences (txt)"))

    args = parser.parse_args()
    exclude_tm_long_gap(args)
