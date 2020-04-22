#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 10.10.2013
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import sys
from trseeker.seqio.fasta_file import sc_iter_fasta
import argparse
from intervaltree import IntervalTree
from collections import defaultdict

from trseeker.tools.jellyfish_tools import query_kmers, Kmer2tfAPI
jf_api = True
try:
    from jellyfish import jellyfish
except:
    print("Failed: from jellyfish import jellyfish")
    try:
        import jellyfish
    except:
        print("Failed: import jellyfish")
        jf_api = False



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Print cov for fasta')
    parser.add_argument('-i','--fasta', help='Fasta file', required=True)
    parser.add_argument('-j','--jf', help='JF file', required=True)
    parser.add_argument('-o','--output', help='Ouput file', required=True)
    parser.add_argument('-k','--k', help='k', required=True)
    parser.add_argument('-a','--annotation', help='Annotation file', required=False, default=False)
    parser.add_argument('-m','--multu', help='Save as multifasta', required=False, default=False)
    parser.add_argument('-n','--new', help='New jf', required=False, default=True)
    parser.add_argument('-c','--cutoff', help='Minimal coverage to report', required=True)
    parser.add_argument('-s','--onlyhits', help='Save only hits', required=False, default=False)
    args = vars(parser.parse_args())

    annotation_file = args["annotation"]
    multi = args["multu"]
    
    if annotation_file:
        tree = IntervalTree()
        with open(annotation_file) as fh:
            data = fh.readlines()
            for line in data:
                line = line.strip().split()
                s,e = int(line[0]), int(line[1])
                if s > e:
                    s,e = e,s
                tree.addi(s-1,e-1, " ".join(line[2:]))

    input_fasta = args["fasta"]

    only_hits = bool(args["onlyhits"])

    jf_path = args["jf"]

    use_new = bool(args["new"])

    k = int(args["k"])
    c = int(args["cutoff"])

    fh = open(args["output"], "w")

    if jf_api and use_new:
        jf_api = jellyfish.QueryMerFile(jf_path)
        kmer2freq = Kmer2tfAPI(jf_api)


    for sid, seq_obj in enumerate(sc_iter_fasta(input_fasta)):
        print seq_obj.header, "Length:", seq_obj.length
        sequence = seq_obj.sequence.upper()
        n = len(sequence)
        kmers = set()
        for i in xrange(n-k+1):
            kmer= sequence[i:i+k]
            tf = kmer2freq[kmer]
            if tf > 0:
                print i, tf
        raw_input("Next item?")


    #         kmers.add()
    #         kmers.add('A' + sequence[i:i+k-1])
    #         kmers.add('C' + sequence[i:i+k-1])
    #         kmers.add('G' + sequence[i:i+k-1])
    #         kmers.add('T' + sequence[i:i+k-1])
    #         kmers.add(sequence[i+1:i+k] + 'A')
    #         kmers.add(sequence[i+1:i+k] + 'C')
    #         kmers.add(sequence[i+1:i+k] + 'T')
    #         kmers.add(sequence[i+1:i+k] + 'G')

    #     print "Kmers:", len(kmers)
    #     print "Load kmers"
        
    #     kmers = list(kmers)
    #     # step = 1000
    #     # for i in xrange(0, len(kmers), step):
    #     #     kmer2freq.update(query_kmers(jf_path, kmers[i:i+step], both_strands=True, verbose=False))
    #     #     print "Loaded", i, i+step, "total size:", len(kmer2freq)

    #     if not (jf_api and use_new):
    #         kmer2freq = defaultdict(int)
    #         kmer2freq = query_kmers(jf_path, kmers, both_strands=True, verbose=30, new=use_new, batch_size=2000)

    #     print "Results:"
    #     for i in xrange(n-k+1):
    #         kmer = sequence[i:i+k].upper()
    #         a_kmer = kmer2freq['A' + kmer[:-1]]
    #         c_kmer = kmer2freq['C' + kmer[:-1]]
    #         g_kmer = kmer2freq['G' + kmer[:-1]]
    #         t_kmer = kmer2freq['T' + kmer[:-1]]
            
    #         kmer_a = kmer2freq[kmer[1:] + 'A']
    #         kmer_c = kmer2freq[kmer[1:] + 'C']
    #         kmer_g = kmer2freq[kmer[1:] + 'G']
    #         kmer_t = kmer2freq[kmer[1:] + 'T']

    #         if kmer_a < c: kmer_a = 0
    #         if kmer_c < c: kmer_c = 0
    #         if kmer_g < c: kmer_g = 0
    #         if kmer_t < c: kmer_t = 0

    #         if a_kmer < c: a_kmer = 0
    #         if c_kmer < c: c_kmer = 0
    #         if g_kmer < c: g_kmer = 0
    #         if t_kmer < c: t_kmer = 0

    #         prev_hits = len([x for x in [a_kmer, c_kmer, g_kmer, t_kmer] if x >= c])
    #         next_hits = len([x for x in [kmer_a, kmer_c, kmer_g, kmer_t] if x >= c])

    #         prefix = "  "
    #         if prev_hits > 1 and next_hits > 1:
    #             prefix = "FF"
    #         elif prev_hits > 1:
    #             prefix = "F-"
    #         elif next_hits > 1:
    #             prefix = "-F"
                
    #         if int(kmer2freq[kmer]) == 0:
    #             if only_hits:
    #                 continue

    #             prefix = "--"
    #             s = prefix, i, "-", kmer, "Prev:", a_kmer, c_kmer, g_kmer, t_kmer, "Next:", kmer_a, kmer_c, kmer_g, kmer_t, "*"*prev_hits, "*"*next_hits, " ".join([x.data for x in tree[i]])
    #         else:
    #             s = prefix, i, kmer2freq[kmer], kmer, "Prev:", a_kmer, c_kmer, g_kmer, t_kmer, "Next:", kmer_a, kmer_c, kmer_g, kmer_t, "*"*prev_hits, "*"*next_hits, " ".join([x.data for x in tree[i]])

    #         print s
    #         fh.write("%s\n" % "\t".join(map(str, s)))

    #     if not multi:
    #         break
    #     

    # fh.close()



