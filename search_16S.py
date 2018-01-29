"""
search_16S.py - a python implementation of the 16S identification method described by
Robert Edgar. See preprint: https://www.biorxiv.org/content/biorxiv/early/2017/04/04/124131.full.pdf

Usage:
    search_16S.py train -i <ref_seqs> -d <ref_db> -k=<kmer_size>
    search_16S.py find -d <ref_db> -i <sample_seqs> [-o <output>]
"""
import skbio
import collections
import docopt
import json

motif_F = "GNTTGATCNTGNC"
motif_R = "AGTCNNAACAAGGTANCNNTA"
max_n_mismtach_in_motif = 4
length_range_for_reporting = (1000, 2500)

def train(ref_seqs, kmer_size, ref_db):
    seqs = skbio.io.read(ref_seqs, format='fasta')
    kmer_counts = collections.Counter()

    for seq in seqs:
        kmer_counts.update(seq.kmer_frequencies(kmer_size))

    total_count = sum(kmer_counts.values())
    kmer_counts = {k: v/total_count for k, v in kmer_counts.items()}
    with open(ref_db, 'w') as fout:
        json.dump(kmer_counts, fout)

def find(sample_seqs, ref_db, outfile):
    with open(ref_db) as fin:
        kmer_counts = json.load(fin)
    kmer_size = len(next(iter(kmer_counts)))
    chance_of_any_kmer = len(kmer_counts) / 4 ** kmer_size

    seqs_to_assess = skbio.io.read(sample_seqs, format='fasta')

    for seq in seqs_to_assess:
        for window in seq.

def main(args):
    if args["train"]:
        train(args["ref_seqs"], args["kmer_size"], args["ref_db"])
    elif args["find"]:
        find(args["sample_seqs"], args["ref_db"], args["output"])

if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)
