#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 19:42:12 2023

@author: sai
"""

import argparse
import time
from multiprocessing.pool import Pool
import os
import sys
import re

import helper
from BloomFilter import BloomFilter

# Parser
parser = argparse.ArgumentParser(description='Microbial Query processor')
parser.add_argument('-m', '--minimizer',
                    help='Minimizer (W:Winnow/T:Topk)', default='W', type=str)
parser.add_argument('-d', '--data_dir',
                    help='data directory', default='/home/sai/summer_school/data/', type=str)
parser.add_argument('-s', '--sim_threshold',
                    help='similarity threshold', default=0.7, type=float)
parser.add_argument('-n', '--max_reads',
                    help='maximum number of reads/species (-1: all reads) ', default=-1, type=int)
parser.add_argument('-i', '--index_file',
                    help='index file', default='index_21mer_640_wm-10_3.bin', type=str)
args = parser.parse_args()

QUERY_FASTQ_DIR = 'query/'
FASTQ_FILE_PREFIX = '.fq'
INDEX_DIR = args.data_dir+'index/'
RESULTS_DIR = args.data_dir+'results/'+QUERY_FASTQ_DIR
INDEX_FILE = INDEX_DIR+args.index_file
TOPK = 5

def run_query_batch(species_id: list) -> dict:
    """
    Receives a list of query species ids and processes each species concurrently by computing:
        a) set of minimizers for all reads in the species
        b) sketch for each read in the species
    and c) Top-k species in the index.

    Returns dictionary containing Top-k species for each read under each species.

    Parameters
    ----------
    species_id : list of query species ids

    Returns
    -------
    dict : (key - species id, values : list containing Top-k tuples.
            Tuple contains similarity score and corresponding species id

    """
    batch_results = {}
    start_time = time.time()
    # Step 1. Read the query fasta files
    if args.minimizer == 'W':
        query_minimizers, read_ids = helper.query_winnow_minimizers(
            f'{args.data_dir+QUERY_FASTQ_DIR}{species_id}{FASTQ_FILE_PREFIX}',
            KMER, THRESHOLD, max_reads=args.max_reads)
    else:
        query_minimizers, read_ids = helper.query_topk_minimizers(
            f'{args.data_dir+QUERY_FASTQ_DIR}{species_id}{FASTQ_FILE_PREFIX}',
            KMER, THRESHOLD, max_reads=args.max_reads)

    print(f"{os.getpid()} : {species_id} {len(query_minimizers)} \
          Minimizers computed in {time.time()-start_time} secs ..")

    # Step 2. Build Species level sketch for each of the read using bloom-filters
    read_results = {}
    log_cnt = 200
    for i in range(len(read_ids)):
        read_id = read_ids[i]
        query_read_minimizer = query_minimizers[i]

        read_sketch = helper.compute_sketch(
            query_read_minimizer, SKETCH_SIZE, NUM_HASHES)

        # query the whole index from ROOT to leaf nodes with threshold
        matches = query_iter(0, read_sketch)
        sorted_matches = sorted(matches, reverse=True)

        # Store only Topk species
        read_results[read_id] = sorted_matches[:TOPK]

        if(i > 0 and i % log_cnt == 0):
            print(
                f' Processed {i} reads for {species_id} in {time.time()-start_time} secs ..')

    batch_results[species_id] = read_results
    return batch_results


def calculate_sim_score(query_sketch: BloomFilter, index_sketch: BloomFilter) -> float:
    """
    Receives Query and Index Sketches and computes match/similarity score between them

    Parameters
    ----------
    query_sketch : BloomFilter
    index_sketch : BloomFilter

    Returns
    -------
    float : match/similarity score.

    """
    match_sketch = query_sketch.intersection(index_sketch)
    match_score = match_sketch.count()[0]/query_sketch.count()[0]
    return match_score


def query_iter(root: int, query_sketch: BloomFilter) -> list:
    """
    Computes the similarity score of the query by traversing the index
    from ROOT nodes to the leaf nodes above a defined threshold iteratively.

    Returns the list of nodes that match the query in the index above a defined threshold.

    Parameters
    ----------
    root : int
    query_sketch : BloomFilter

    Returns
    -------
    list

    """
    node_idxs = []
    node_idxs.append(root)
    matches=[]
    while len(node_idxs) > 0:
        i = node_idxs.pop()
        if node_idxs == -1:
            continue
        lnode_idx, index_sketch = sketch_index[i]
        sim_score = calculate_sim_score(query_sketch, index_sketch)
        if sim_score >= args.sim_threshold:
            if lnode_idx != -1:
                node_idxs.append(lnode_idx)
                node_idxs.append(lnode_idx-1)
            else:
                matches.append((sim_score, IDX_ID[i]))
    return matches


def compute_metrics(query_results: dict) -> (list, list):
    """
    Computes metrics viz., Recall at defined Top-k across all the species.
    Returns a) Recall@Top-k as list ,and b) True positives @Top-k

    Parameters
    ----------
    query_results : dict
        DESCRIPTION.

    Returns
    -------
    list : Recall@Top-k
    list : True positives @Top-k
    """
    true_positives = [0]*TOPK
    nreads=0
    for spid in query_results:
        read_results = query_results[spid]
        nreads += len(read_results)
        for read_id in read_results:
            topk_species = read_results[read_id]
            topk_len=len(topk_species)
            for i in range(TOPK):
                true_positives[i] += 1 if(topk_len-1>=i and topk_species[i][1] == spid) else 0

    #print(true_positives)
    recall = [0.0]*TOPK
    for i in range(TOPK):
        if i > 0: #cummulative
            true_positives[i] += true_positives[i-1]
        recall[i] = round(true_positives[i]/nreads, 4)
    return recall, true_positives

if __name__ == '__main__':
    global IDX_ID, ID_IDX, sketch_index, SKETCH_SIZE, NUM_HASHES, KMER, THRESHOLD

    t0 = time.time()
    FILE_PATTERN = '(index_)([0-9]+)(mer_)([0-9]+)_(wm|tm)-([0-9]+)_([0-9]+).bin'

    if not os.path.isfile(INDEX_FILE):
        print(f"ERROR : {INDEX_FILE} doesn't exist! Please check the path and run again ..")
        sys.exit()

    pattern_m = re.match(FILE_PATTERN, args.index_file)
    if pattern_m is None :
        print(
            f'ERROR : {INDEX_FILE} invalid! Please recreate the index file using index.py')
        sys.exit()

    KMER, SKETCH_SIZE, THRESHOLD, NUM_HASHES = int(pattern_m.group(2)), int(
        pattern_m.group(4)), int(pattern_m.group(6)), int(pattern_m.group(7))

    print(
        f' Minimizer : {args.minimizer}, Kmer : {KMER}, Sketch size : {SKETCH_SIZE}')
    print(f' Winnow Threshold : {THRESHOLD}, Hashes : {NUM_HASHES}')
    print(f' FASTQ DIR : {QUERY_FASTQ_DIR}, RESULTS DIR : {RESULTS_DIR} ')
    print(args)

    # Step 1. Load the index from file
    sketch_index = helper.read_index(
        INDEX_FILE, sketch_size=SKETCH_SIZE, num_hashes=NUM_HASHES)
    assert len(sketch_index) == 2047
    print(f' Loading {len(sketch_index)} Index took {time.time()-t0} secs')

    # Step 2. get index to species id mapping
    TREE_FILE = INDEX_DIR + \
        args.index_file.replace('bin', 'csv').replace('index', 'tree')
    IDX_ID, ID_IDX = helper.get_idx_to_id_map(TREE_FILE)

    # Step 3. Read species fasta files from the data directory
    query_species_ids = [line.replace(FASTQ_FILE_PREFIX, '')
                         for line in os.listdir(args.data_dir+QUERY_FASTQ_DIR)]
    assert len(query_species_ids) > 0
    print(f" Num. Query Species : {len(query_species_ids)}")

    print(' Building Bloom Filters .. ')
    # Step 3. compute the sketches for each of the query species
    query_results = {}
    with Pool() as pool:
        #bt = time.time()
        for batch_results in pool.imap(run_query_batch, query_species_ids):
            query_results.update(batch_results)
            # print(
            #     f" Processed Batch of {len(batch_results)} in {time.time()-bt} secs ..")

    # More than 1 species is required to compute metrics
    assert len(query_results) > 1
    recall_topk, _ = compute_metrics(query_results)
    if not os.path.isdir(RESULTS_DIR):
        os.mkdir(RESULTS_DIR)

    RESULTS_FILE = RESULTS_DIR + \
        args.index_file.replace('bin', 'csv').replace('index', 'results')
    with open(RESULTS_FILE, 'w', encoding='utf-8') as resf:
        resf.write("\t".join('Rec@'+str(i+1) for i in range(len(recall_topk)))+'\n')
        for i in range(len(recall_topk)):
            resf.write(str(recall_topk[i])+'\t')
        resf.write('\n')

    print(recall_topk)
    print(f' Took {time.time()-t0} seconds to complete the pipeline..')
    
