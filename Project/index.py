#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 11:08:29 2023

@author: sai
"""
import argparse
import time
from multiprocessing.pool import Pool
import os

import helper
from Node import Node

# Parser
parser = argparse.ArgumentParser(description='Microbial Indexer')
parser.add_argument('-m', '--minimizer',
                    help='Minimizer (W:Winnow/T:Topk)', default='W', type=str)
parser.add_argument('-t', '--topk_threshold',
                    help='Topk threshold Minimizer', default=2000, type=int)
parser.add_argument('-w', '--winnow_threshold',
                    help='Winnow threshold Minimizer', default=50, type=int)
parser.add_argument('-k', '--kmer',
                    help='kmer length', default=21, type=int)
parser.add_argument('-s', '--sketch_size',
                    help='sketch size in KiloBytes', default=640, type=int)
parser.add_argument('-n', '--num_hashes',
                    help='sketch size in KiloBytes', default=3, type=int)
parser.add_argument('-d', '--data_dir',
                    help='data directory', default='/home/sai/summer_school/data/', type=str)

args = parser.parse_args()

FASTA_DIR = 'input/'
FASTA_FILE_PREFIX = '.fa'
PARAM_STR = f'{args.kmer}mer_{args.sketch_size}_{"wm-"+str(args.winnow_threshold) if(args.minimizer=="W") else "tm-"+str(args.topk_threshold)}_{args.num_hashes}'
INDEX_DIR = args.data_dir+'index/'
INDEX_FILE = f'index_{PARAM_STR}.bin'
TREE_FILE = f'tree_{PARAM_STR}.csv'

def run_sketch_batch(species_id: list) -> dict:
    """
    

    Parameters
    ----------
    species_id : list
        DESCRIPTION.

    Returns
    -------
    dict
        DESCRIPTION.

    """
    st = time.time()
    batch_sketches = {}
    # Step 1. Read the Fastq file and Build Minimizers

    if args.minimizer == 'W':
        # If winnow threshold is set to greater than 0, use winnowing minimizer
        minimizers = helper.compute_winnowing_minimizers(
            f'{args.data_dir+FASTA_DIR}{species_id}{FASTA_FILE_PREFIX}',
            args.kmer, args.winnow_threshold)
    else:
        # use Topk-minimizer
        minimizers = helper.compute_topk_minimizers(
            f'{args.data_dir+FASTA_DIR}{species_id}{FASTA_FILE_PREFIX}',
            args.kmer, args.topk_threshold)

    #print(f"pid: {os.getpid()}, species_id : {species_id} minimizers : {len(minimizers)} \
    #    Minimizers computed in {time.time()-st} secs ..")
    #helper.writeMinimizers(minimizers, "index_minimizers")

    # Step 2. Build Species level sketch using bloom-filters
    sketch = helper.compute_sketch(
        minimizers, args.sketch_size, args.num_hashes)
    batch_sketches[species_id] = sketch

    return batch_sketches

def run_sketch(sp_ids: list) -> dict:
    """
    

    Parameters
    ----------
    sp_ids : list
        DESCRIPTION.

    Returns
    -------
    dict
        DESCRIPTION.

    """
    sk_dict = {}
    with Pool() as pool:
        #bt = time.time()
        for batch_sketches in pool.imap(run_sketch_batch, sp_ids):
            sk_dict.update(batch_sketches)
            #print(f" Processed Batch in {time.time()-bt} secs ..")
    
    return sk_dict

def build_tree(sk_dict: dict) -> list:
    """
    

    Parameters
    ----------
    sk_dict : dict
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    """
    sk_nodes = []
    # create sketches as Node in the alphabetical order
    sorted_species_id = sorted(sk_dict.keys())
    i = 0
    node_len = len(sorted_species_id)
    # Build node_idx from bottom up (leaf nodes-root node)
    node_idx = (2*node_len)-2
    while i < node_len:
        sk_nodes.append(
            Node(sk_dict[sorted_species_id[i]], node_idx, sorted_species_id[i]))
        i += 1
        node_idx -= 1

    print(f' Initial Node sketches : {len(sk_nodes)} ')

    # construct a perfect binary tree from leaf to root
    sketch_q = sk_nodes.copy()
    while len(sketch_q) > 1:
        lnode = sketch_q.pop(0)
        rnode = sketch_q.pop(0)
        # compute the union of left and right node sketches
        merged_node_val = lnode.val.union(rnode.val)
        #print(f' merged node sketch count : {merged_node_val.count()}')

        # create a new merged node and add it to the queue
        sketch_q.append(Node(merged_node_val, node_idx,
                        'IN', left=lnode, right=rnode))
        node_idx -= 1

        # Add the new node to the final list
        sk_nodes.append(sketch_q[-1])

    return sk_nodes


if __name__ == '__main__':

    t0 = time.time()

    print(args)

    # Read species fasta files from the data directory
    species_ids = [line.replace(FASTA_FILE_PREFIX, '')
                   for line in os.listdir(args.data_dir+FASTA_DIR)]
    print(f" Total Species : {len(species_ids)}")

    # Step 1. Compute sketches for each species in the dataset
    sketch_dict = run_sketch(species_ids)
    print(f' Processed {len(sketch_dict)} Species ...')

    # Step 2. Construct Tree using the sketches
    sketch_nodes = build_tree(sketch_dict)
    print(f' Final Node sketches : {len(sketch_nodes)} ')
    helper.write_sketch_tree(sketch_nodes, INDEX_DIR+TREE_FILE)

    # Step 3. Write the sketches (from leaf nodes to ROOT node) to file
    helper.write_index(sketch_nodes, INDEX_DIR+INDEX_FILE)

    print(f' Took {time.time()-t0} seconds to complete the pipeline..')
