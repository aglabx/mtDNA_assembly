#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 10.10.2013
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import sys
from trseeker.seqio.fasta_file import sc_iter_fasta
import argparse

from trseeker.tools.jellyfish_tools import query_kmers, Kmer2tfAPI
from trseeker.tools.ngrams_tools import print_next_cutoff
try:
    from Ariadna.AriadnaPy.tools.aindex import load_aindex, get_left_right_distances
except:
    load_aindex = None
    get_left_right_distances = lambda x: x

from collections import defaultdict

from trseeker.tools.sequence_patterns import DNA2IUPAC

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


def extend_sequence_right(kmer, kmer2freq, length=23):
    
    sequence = kmer
    global_i = 0
    while True:
        kmer_a = kmer[1:] + 'A'
        kmer_c = kmer[1:] + 'C'
        kmer_g = kmer[1:] + 'G'
        kmer_t = kmer[1:] + 'T'
        R, n, nucleotides, max_hits = print_next_cutoff(kmer, kmer2freq, cutoff=c)
        
        if n == 1:
            pos = -1
        else:
            return None
        sequence += max_hits[pos][1]
        kmer = sequence[-k:]
        global_i += 1;

        if global_i == length:
            return sequence[k-1:]


# def check_snp_bubble(kmer1, kmer2, kmer2freq, length=23):

#     print kmer1, extend_sequence_right(kmer1, kmer2freq, length, k=k)
#     print kmer2, extend_sequence_right(kmer2, kmer2freq, length, k=k)

def check_pairs(left_kmer, kmer_a, kmer_c, kmer_t, kmer_g, aindex, sequence, key_kmers_solved, nucleotides, key_kmers_unsolved, shift=None, k=23):
    count_a = len(get_left_right_distances(left_kmer,kmer_a,aindex))
    count_c = len(get_left_right_distances(left_kmer,kmer_c,aindex))
    count_t = len(get_left_right_distances(left_kmer,kmer_t,aindex))
    count_g = len(get_left_right_distances(left_kmer,kmer_g,aindex))

    print "Shift:", shift, max_hits
    print "A", count_a, extend_sequence_right(kmer_a, kmer2freq, length=k)
    print "C", count_c, extend_sequence_right(kmer_c, kmer2freq, length=k)
    print "T", count_t, extend_sequence_right(kmer_t, kmer2freq, length=k)
    print "G", count_g, extend_sequence_right(kmer_g, kmer2freq, length=k)

    if count_a > min_error and count_c < min_error and count_t < min_error and count_g < min_error:
        sequence += 'A'
        key_kmers_solved[kmer_a] = len(sequence) - k
        for nucl in nucleotides:
            if nucl == 'A':
                continue
            key_kmers_unsolved[kmer_a].append(kmer_a[:-1]+nucl)
        return sequence
    if count_c > min_error and count_a < min_error and count_t < min_error and count_g < min_error:
        sequence += 'C'
        key_kmers_solved[kmer_c] = len(sequence) - k
        for nucl in nucleotides:
            if nucl == 'C':
                continue
            key_kmers_unsolved[kmer_c].append(kmer_c[:-1]+nucl)
        return sequence
    if count_t > min_error and count_c < min_error and count_a < min_error and count_g < min_error:
        sequence += 'T'
        key_kmers_solved[kmer_t] = len(sequence) - k
        for nucl in nucleotides:
            if nucl == 'T':
                continue
            key_kmers_unsolved[kmer_t].append(kmer_t[:-1]+nucl)
        return sequence
    if count_g > min_error and count_c < min_error and count_t < min_error and count_a < min_error:
        sequence += 'G'
        key_kmers_solved[kmer_g] = len(sequence) - k
        for nucl in nucleotides:
            if nucl == 'G':
                continue
            key_kmers_unsolved[kmer_g].append(kmer_g[:-1]+nucl)
        return sequence

    return None

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Print cov for fasta')
    parser.add_argument('-i','--sequence', help='Start sequence', required=True)
    parser.add_argument('--right', help='Stop sequence', required=False, default=None)
    parser.add_argument('-j','--jf', help='JF file', required=True)
    parser.add_argument('-k','--k', help='k', required=False, default=23)
    parser.add_argument('-n','--new', help='Use new', required=False, default=True)
    parser.add_argument('-s','--start', help='start', required=False, default=0)
    parser.add_argument('--min_error', help='error cutoff', required=False, default=0)
    parser.add_argument('-p','--pause', help='pause', required=False, default=False)
    parser.add_argument('-x','--xkmer', help='Remove seed kmer', required=False, default=False)
    parser.add_argument('-o','--output', help='output', required=False, default=False)
    parser.add_argument('--index', help='index', required=False, default=False)
    parser.add_argument('--reads', help='reads', required=False, default=False)
    parser.add_argument('--aindex', help='aindex', required=False, default=False)
    parser.add_argument('-c','--cutoff', help='Minimal coverage to report', required=False, default=0)
    args = vars(parser.parse_args())

    input_fasta = args["sequence"]

    stop_kmer = args["right"]

    index_prefix = args["index"]
    reads_file = args["reads"]
    aindex_prefix = args["aindex"]

    xlengths = len(input_fasta)

    start = int(args["start"])

    jf_path = args["jf"]

    k = int(args["k"])
    c = int(args["cutoff"])
    pause = args["pause"]
    min_error = int(args["min_error"])

    xkmer = args["xkmer"]

    use_new = bool(args["new"])

    output_file = args["output"]
    
    sequence = input_fasta.upper()
    kmer = sequence[-k:]

    ref_kmer = kmer

    global_i = 0


    if min_error:
        settings = {
            "index_prefix": index_prefix,
            "aindex_prefix": aindex_prefix,
            "reads_file": reads_file,
            "max_tf": 10000000,
        }
        aindex = load_aindex(settings)

    if output_file:
        fh = open(output_file, "w")

    if output_file:
        fh.write(sequence)

    
    if jf_api and use_new and jf_path.endswith(".jf2"):
        jf_api = jellyfish.QueryMerFile(jf_path)
        kmer2freq = Kmer2tfAPI(jf_api)
    else:
        import aindex
        settings = {
            "index_prefix": jf_path,
        }
        kmer2freq = aindex.load_aindex(settings, skip_aindex=True, skip_reads=True)



    key_kmers_solved = {}
    key_kmers_unsolved = defaultdict(list)

    while True:

        if stop_kmer and sequence.endswith(stop_kmer):
            print "Solution found"
            print sequence
            break


        kmer_a = kmer[1:] + 'A'
        kmer_c = kmer[1:] + 'C'
        kmer_g = kmer[1:] + 'G'
        kmer_t = kmer[1:] + 'T'

        if global_i and ref_kmer == kmer:
            print "Found start kmer"
            break

        if not use_new:
            kmer2freq = query_kmers(jf_path, [kmer_a, kmer_t, kmer_g, kmer_c], both_strands=True, verbose=False)

        R, n, nucleotides, max_hits = print_next_cutoff(kmer, kmer2freq, cutoff=c)

        # print kmer
        # print R, n, nucleotides, max_hits
        # raw_input("?")

        print n, max_hits, nucleotides

        global_i += 1;
            
        if n == 1:
            if pause:
                raw_input("Next?")
            pos = -1

            sequence += max_hits[pos][1]

        else:
            shift = 0
            full_stop = False 
            while True:
                if min_error:

                    left_kmer = sequence[-40-shift:-40+k-shift]
                    result = check_pairs(left_kmer, kmer_a, kmer_c, kmer_t, kmer_g, aindex, sequence, key_kmers_solved, nucleotides, key_kmers_unsolved, shift=None, k=k)

                    if result is not None:
                        sequence = result
                        break

                pos = raw_input("Fork (or STOP)?") or None
                # with open("right.temp", "w") as fht:
                #     fht.write(sequence)
                


                if pos.strip() == "STOP":
                    full_stop = True
                    break

                if pos != None:
                    pos = pos.strip()
                    try:
                        pos = int(pos)
                        if pos > 4:
                            shift = pos
                            continue
                        sequence += max_hits[pos][1]
                        key_kmers_solved[sequence[-k:]] = len(sequence) - k
                    except:
                        if len(pos) > 1:
                            if pos in DNA2IUPAC:
                                sequence += DNA2IUPAC[pos] + extend_sequence_right(kmer_a[:-1]+pos[0], kmer2freq, length=k)[1:]
                                break
                            else:
                                sequence += pos
                                key_kmers_solved[sequence[-k:]] = len(sequence) - k
                                global_i += len(pos)-1
                                break    
                        if pos and pos.strip() in ["A","C","T","G"]:
                            sequence += pos
                            key_kmers_solved[sequence[-k:]] = len(sequence) - k
                            break
                        continue
                    break

                # if min_error:
                #     print "Check solved:"
                #     keys = key_kmers_solved.keys()
                #     while keys:
                #         left_kmer = keys.pop()
                #         print "Positive:"
                #         pos_result  = check_pairs(left_kmer, kmer_a, kmer_c, kmer_t, kmer_g, aindex, sequence, key_kmers_solved, nucleotides, key_kmers_unsolved, shift=None, k=23)
                #         print "Negative:"
                #         for neg_kmer in key_kmers_unsolved[left_kmer]:
                #             neg_result = check_pairs(left_kmer, kmer_a, kmer_c, kmer_t, kmer_g, aindex, sequence, key_kmers_solved, nucleotides, key_kmers_unsolved, shift=None, k=23)

                #         if result is not None:
                #             raw_input("Done?")
                #             sequence = result
                #             break
                #         raw_input("Miss?")

                shift += 5
        
            if full_stop:
                break

        kmer = sequence[-k:]
        if output_file:
            print start, start+global_i, len(sequence), sequence[-53:], start+global_i+1+k
        else:
            if xkmer:
                print start, start+global_i, len(sequence), sequence[xlengths:]
            else:
                print start, start+global_i, len(sequence), sequence

        if output_file:
            fh.write(sequence[-1])

    if output_file:
        fh.close()
