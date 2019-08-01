#!/usr/bin/env python
import argparse

from Bio import AlignIO


# 3. The list of TMs is sorted by start position. Then j, which represents the
#    the position in the original reference sequence is set to -1.
#    A for loop then goes through aln_seq with two iterators, i and x. i
#    represents the position in the alignment, and x represents the amino acid
#    at position i.
#    j only increases if x is not a gap "-". When j is equal to the start
#    position of of the first TM, i is recorded as the start position of the
#    TM1 in the alignment and returned.
def find_tm1(aln_seq, or_tms):
    or_tms.sort(key=lambda x: x[1])
    j = -1
    for i, x in enumerate(aln_seq):
        if x != "-":
            j += 1
            if j == or_tms[0][1]:
                tm1_start = i
                break
    return tm1_start


# 2. Takes in arguments and makes:
#    ref_or a SeqIO Seqrecord,
#    orf_aln an AlignIO alignment
#    in_ref_or_tm an array of (Ref OR ID, TM start, TM end, -1, -1)
#    rem_id an array
#    Then it checks that all TMs have the same chrom value. If not, error
#    The function then attempts to get the index of the reference or. If it
#    can't, error.
#    The alignment sequence of REF OR, and the list of TM regions are sent to
#    find_tm1.
def assign_start_codon(args):
    with open(args.orf_aln, "rU") as in_orf_aln:
        orf_aln = AlignIO.read(in_orf_aln, "fasta")
    with open(args.rem_id, "rU") as in_rem_id:
        remove_names = in_rem_id.read().splitlines()
    or_tms = []
    ref_id = []
    with open(args.ref_or_tm, "rU") as in_ref_or_tm:
        for line in in_ref_or_tm:
            cols = line.rstrip().split()
            assert(len(cols) == 3)
            # Ref_OR ID, Ref TM Start, Ref TM Stop
            or_tms.append([cols[0], int(cols[1])])
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

    tm1_start = find_tm1(orf_aln[ref_or_index].seq, or_tms)
# 4. The region from the start of the sequence to the start of TM1 is
#    extracted from the alignment as upstream. If upstream contains a start
#    codon, then it divides upstream into 3 regions: a (start of ORF to -34)
#    b (-34 to -20) and c (-20 to TM1). Depending on whether a, b or c
#    contain start codons, the index of the start codon is decided. This index
#    is either the last M in b, the last M in a, or the first m in c. then
#    each orf is formatted as a gff entry with its new start codon.
    upstream = orf_aln[:, :tm1_start]
    real_start = []
    for record in upstream:
        rec_start = 0
        if record.id != ref_id:
            if record.seq.count("M") < 1:
                remove_names.append(record.id)
            elif record.seq.count("M") > 1:
                up_seq = str(record.seq).replace("-", "")
                up_len = len(up_seq)
                a, b, c = "", "", ""
                if up_len > 20:
                    c = up_seq[-20:]
                    if up_len > 34:
                        b = up_seq[-34:-20]
                        a = up_seq[:-34]
                    else:
                        b = up_seq[:-20]
                else:
                    c = up_seq
                if "M" in b:
                    rec_start = len(a) + b.rfind("M")
                elif "M" in a:
                    rec_start = a.rfind("M")
                else:
                    rec_start = len(a) + len(b) + c.find("M")
        real_start.append(rec_start)

    for i in range(0, len(orf_aln)):
        if orf_aln[i].id not in remove_names:
            chrom = orf_aln[i].id.split(":")[0]
            strand = orf_aln[i].id.split("(")[1].split(")")[0]
            ref_start = int(
                orf_aln[i].id.split(":")[1].split("(")[0].split("-")[0]) + 1
            ref_end = int(
                orf_aln[i].id.split(":")[1].split("(")[0].split("-")[1])
            if strand == "+":
                ref_start += real_start[i] * 3
            else:
                ref_end -= real_start[i] * 3
            # output (GFF format)
            orf = [chrom, "tm_gap", "ORF",
                   str(ref_start), str(ref_end), "0.0", strand, ".", "."]
            print('\t'.join(orf))


# 1. Takes in arguments genome, orf_aln, or_ref, or_ref_tm and
#    hands them to assign_start_codon
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assignment of a proper initiation codon")
    parser.add_argument(
        '--orf_aln',
        required=True,
        help=("OR reference sequence (Results of multiple alignment)"))
    parser.add_argument(
        '--ref_or_tm',
        required=True,
        help=("TMs of OR reference sequence (BED)"))
    parser.add_argument(
        '--rem_id',
        required=True,
        help=("ids of reference or sequences (txt)"))

    args = parser.parse_args()
    assign_start_codon(args)
