import argparse

from Bio import SeqIO


<<<<<<< Updated upstream
# 3. If seq is on the negative strand, the reverse complement of seq is
#    assigned to nuc. orffinder then goes through nuc 3 times to get all
#    reading frames on strand.
#    Each potential reading frame is translated
#    If any of the following steps fails the program stops looking in the
#    potential reading frame
#    The location of the first "M" in the translated sequence is found
#    The first stop codon ("*") after that "M" is found
#    The resulting protein is at least min_protein_length AA long
#    There are no unkown proteins ("X") in the sequence
#    the start and stop positions on the original sequence are recorded
def orffinder(seq, min_protein_length, strand):
=======
# 6. orffinder goes through seq 6 times to get all 6 possible reading frames
#    Each time the potential reading frame is translated
#    A while loop goes through each amino acid in the reading frame.
#    If an AA is methionine then the next stop codon "*" is located
#    If the resulting fragment is longer than min protein length, then the
#    start, stop and strand are added to orfs
#    After possible sequences from all 6 reading frames are added, orfs is
#    sorted and returned.
def orffinder(seq, min_protein_length):
>>>>>>> Stashed changes
    orfs = []
    seq_len = len(seq)
    if strand == "-":
        nuc = seq.reverse_complement()
    else:
        nuc = seq
    for frame in [0, 1, 2]:
        # seq_len 300, trans_len: 300, 297, 297
        trans_len = ((seq_len - frame) // 3) * 3
        # nuc[0:300], [1:298], [2:299]
        trans = str(nuc[frame:frame + trans_len].translate()).upper()
        aa_len = len(trans)
        aa_start = 0
        aa_stop = 0
        while aa_start < aa_len:
            # finds first M in nuc aa_start = 24
            aa_start = trans.find("M", aa_start)
            if aa_start == -1:
                break
            # finds first * after M aa_stop = 74
            aa_stop = trans.find("*", aa_start)
            if aa_stop == -1:
                break
            # Checks if protein is long enough
            if aa_stop - aa_start >= min_protein_length:
                # Rejects proteins with unknown AA
                if "X" in trans[aa_start:aa_stop]:
                    break
                if strand == "+":
                    # dna start = 73
                    dna_start = frame + aa_start * 3
                    # dna stop =  226
                    dna_stop = frame + (aa_stop + 1) * 3
                else:
                    #
                    dna_start = seq_len - (frame + (aa_stop + 1) * 3)
                    dna_stop = seq_len - (frame + aa_start * 3)

                orfs.append((dna_start + 1, dna_stop, strand))
            aa_start = aa_stop + 1
    orfs.sort()
    return orfs

<<<<<<< Updated upstream

# 2. Defines  min_protein_length and converts the subject genome
#    to a dictionary
#    For each hit in the best-hit file, extracts the corresponding sequence
#    from the subject genome and sends it and min_protein lenghth to orffinder
# 4. Each potential reading frame is added to orf, which is then printed as a
#    BED file
=======

# 4. Extends seq upstream and downstream flank bp or to the start/end of the
#    reference sequence, whichever is shortest
def getsubseq(seq, start, end, flank):
    seq_len = len(seq)
    seq_start = max(0, start - flank)
    seq_end = min(end + flank, seq_len)
    return seq_start, seq_end, seq[seq_start - 1:seq_end]


# 2. Defines flank and min_protein_length and converts the reference sequence
#    to a dict
# 3. For each hit in the best-hit file, extracts the corresponding reference
#    sequence and sends it, flank and the start and end positions of the hit to
#    getsubseq
# 5. Send extended sequence and min_protein_length to orffinder
# 6. Each potential reading frame is added to orf, which is then printed
>>>>>>> Stashed changes
def findorf(args):
    min_protein_length = args.min_p_len
    with open(args.subject, "rU") as in_subject:
        subj_dict = SeqIO.to_dict(SeqIO.parse(in_subject, "fasta"))
    with open(args.hit, "rU") as in_hit:
        for line in in_hit:
            hit = line.rstrip().split()
            chrom = hit[0]
            subj_record = subj_dict[chrom]
            subj_start, subj_stop = int(hit[3]), int(hit[4])
            strand = hit[6]
            seq = subj_record.seq[subj_start - 1:subj_stop]
            orfs = orffinder(seq, min_protein_length, strand)
            for orf_start, orf_stop, strand in orfs:
                orf = [subj_record.id, "findorf", "ORF",
                       str(subj_start - 1 + orf_start),
                       str(subj_start - 1 + orf_stop),
                       "0.0", strand, ".", "."]
                print('\t'.join(orf))


# 1. Defines required arguments and calls findorf
if __name__ == '__main__':
    parser = argparse.ArgumentParser("find open reading frame")
    parser.add_argument(
        "-S", "--subject",
        required=True,
        help=("subject genome (FASTA)"))
    parser.add_argument(
        "-H", "--hit",
        required=True,
        help=("blast hits (GFF)"))
    parser.add_argument(
        "-m", "--min_p_len",
        type=int,
        default=250,
        help=("minimum length of protein [250] (int)"))

    args = parser.parse_args()
    findorf(args)
