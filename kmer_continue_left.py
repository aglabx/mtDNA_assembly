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

from trseeker.tools.jellyfish_tools import query_kmers, Kmer2tfAPI
from trseeker.tools.ngrams_tools import print_prev_cutoff
from trseeker.tools.sequence_tools import get_revcomp

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
    parser.add_argument('-i','--sequence', help='Start sequence', required=True)
    parser.add_argument('-j','--jf', help='JF file', required=True)
    parser.add_argument('-n','--new', help='Use new', required=False, default=True)
    parser.add_argument('-k','--k', help='k', required=False, default=23)
    parser.add_argument('-s','--start', help='start', required=False, default=0)
    parser.add_argument('-c','--cutoff', help='Minimal coverage to report', required=False, default=0)
    parser.add_argument('--upper', help='Maximal     coverage to report', required=False, default=0)
    parser.add_argument('-o','--output', help='output', required=False, default=False)
    args = vars(parser.parse_args())

    input_fasta = args["sequence"]

    start = int(args["start"])

    jf_path = args["jf"]

    k = int(args["k"])
    c = int(args["cutoff"])

    use_new = bool(args["new"])
    sequence = input_fasta.upper()
    kmer = sequence[:k]

    output_file = args["output"]

    print sequence
    print "Kmer:", kmer

    

    if jf_api and use_new:
        jf_api = jellyfish.QueryMerFile(jf_path)
        kmer2freq = Kmer2tfAPI(jf_api)


    global_i = 0

    ref_kmer = kmer

    final_sequence = [kmer]

    while True:

        a_kmer = 'A' + kmer[:-1]
        c_kmer = 'C' + kmer[:-1]
        g_kmer = 'G' + kmer[:-1]
        t_kmer = 'T' + kmer[:-1]

        if not use_new:
            kmer2freq = query_kmers(jf_path, [a_kmer, c_kmer, t_kmer, g_kmer], both_strands=True, verbose=False)

        R, n, nucleotides, max_hits = print_prev_cutoff(kmer, kmer2freq, cutoff=c)

        if global_i and ref_kmer == kmer:
            print "Found start kmer"
            break

        print n, max_hits

        global_i -= 1;
            
        skip = False
        if n == 1:
            # raw_input("Next?")
            pos = -1
        else:
            while True:
                pos = raw_input("Fork (number, save, dec, inc)?") or None
                # with open("left.temp", "w") as fh:
                #     fh.write(sequence)
                if pos.strip() in ["dec", "inc"]:
                    pos = raw_input("Set c:") or 0
                    c = 0
                    if pos != 0:
                        c = int(pos.strip()) 
                    skip = True
                    break
                if pos.strip() in ["save"]:
                    with open(output_file, "w") as fh:
                        fh.write("".join(final_sequence))
                    print "Saved."

                if pos != None:
                    try:
                        pos = int(pos.strip())
                    except:
                        continue
                    break

                
        if skip:
            continue
        sequence = max_hits[pos][1] + sequence
        kmer = sequence[:k]
        if output_file:
            print start, start+global_i, len(sequence), sequence[:53]
        else:
            print start, start+global_i, len(sequence), sequence


        if output_file:
            final_sequence.insert(0, max_hits[pos][1])

    if output_file:
        with open(output_file, "w") as fh:
            fh.write("".join(final_sequence))
        
