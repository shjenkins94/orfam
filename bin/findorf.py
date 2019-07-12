import argparse

from Bio import SeqIO


# 6. orffinder goes through seq 6 times to get all 6 possible reading frames
#    Each time the potential reading frame is translated
#    A while loop goes through each amino acid in the reading frame.
#    If an AA is methionine then the next stop codon "*" is located
#    If the resulting fragment is longer than min protein length, then the 
#    start, stop and strand are added to orfs
#    After possible sequences from all 6 reading frames are added, orfs is 
#    sorted and returned.
def orffinder(seq, min_protein_length):
    orfs = []
    seq_len = len(seq)

    for strand, nuc in [("+", seq), ("-", seq.reverse_complement())]:
        for frame in [0, 1, 2]:
            trans_len = ((seq_len - frame) // 3) * 3
            trans = str(nuc[frame:frame + trans_len].translate()).upper()
            aa_len = len(trans)
            aa_start = 0
            aa_stop = 0
            while aa_start < aa_len:
                aa_start = trans.find("M", aa_start)
                if aa_start == -1:
                    break
                aa_stop = trans.find("*", aa_start)
                if aa_stop == -1:
                    break
                if aa_stop - aa_start >= min_protein_length:
                    if "X" in trans[aa_start:aa_stop]:
                        break
                    if strand == "+":
                        dna_start = frame + aa_start * 3
                        dna_stop = frame + (aa_stop + 1) * 3
                    else:
                        dna_start = seq_len - (frame + (aa_stop + 1) * 3)
                        dna_stop = seq_len - (frame + aa_start * 3)

                    orfs.append((dna_start + 1, dna_stop, strand))

                aa_start = aa_stop + 1
    orfs.sort()
    return orfs

# 4. Extends seq upstream and downstream flank bp or to the start/end of the reference 
#    sequence, whichever is shortest
def getsubseq(seq, start, end, flank):
    seq_len = len(seq)
    seq_start = max(0, start - flank)
    seq_end = min(end + flank, seq_len)
    return seq_start, seq_end, seq[seq_start - 1:seq_end]


# 2. Defines flank and min_protein_length and converts the reference sequence to a dict
# 3. For each hit in the best-hit file, extracts the corresponding reference sequence
#    and sends it, flank and the start and end positions of the hit to getsubseq  
# 5. Send extended sequence and min_protein_length to orffinder
# 6. Each potential reading frame is added to orf, which is then printed
def findorf(args):
    flank = 1000
    min_protein_length = 250
    ref_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    for line in args.hit:
        hit = line.rstrip().split()
        chrom = hit[0]
        ref_record = ref_dict[chrom]
        ref_start, ref_stop, seq = getsubseq(
            ref_record.seq, int(hit[3]), int(hit[4]), flank)
        orfs = orffinder(seq, min_protein_length)
        for orf_start, orf_stop, strand in orfs:
            orf = [ref_record.id, "findorf", "ORF",
                   str(ref_start + (orf_start - 1)),
                   str(ref_start + (orf_stop - 1)),
                   "0.0", strand, ".", "."]
            print('\t'.join(orf))


# 1. Defines required arguments and calls findorf
if __name__ == '__main__':
    parser = argparse.ArgumentParser("find open reading frame")
    parser.add_argument(
        "reference", type=argparse.FileType("r"),
        help=("FASTA    reference genome, FASTA format"))
    parser.add_argument(
        "hit", type=argparse.FileType("r"), help=("GFF  blast output"))

    args = parser.parse_args()
    findorf(args)
