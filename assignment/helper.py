#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 23:32:26 2023

@author: sai
"""
import heapq
import mmh3
from BloomFilter import BloomFilter
from bitarray import bitarray

SEED=42
def check_non_canonicalbases(seq: str) -> (list, int):
    """
    Check the input sequence if there are any Non-canonical bases (below)
    in the sequences.
    Example : BKNDYWMVRH (Non-canonical bases) - NYR (89% - 3969 out of 4496)
    Returns the positions of the non-canonical bases and the number of non-canonical bases
    in the sequence.
    Parameters
    ----------
    seq : str

    Returns
    -------
    (list, int)

    """
    pos = []
    for j in range(len(seq)):
        if(seq[j] != 'G' and seq[j] != 'C' and seq[j] != 'T' and seq[j] != 'A'):
            pos.append(j)
    return pos, len(pos)


def query_winnow_minimizers(fastq_file: str, k: int, winnow_t: int, max_reads=-1) -> set:
    """
    Computes the minimizers based on defined winnow threshold for each read in the query fastq file.
    If max_reads is set to -1, then minimizers are computed for all reads,
    else minimizers are computed upto 'max_reads' reads.

    Query Fastq file format:
    (line1) @read information
    (line2) read sequence
    (line3 and line4) read quality information (ignored)

    Returns Winnow-based minimizers for every query read

    Parameters
    ----------
    fastq_file (query read file) : str
    k : int
    threshold : int
    max_reads : int, optional

    Returns
    -------
    list

    """
    minimizers = []
    read_ids=[]
    # Open the fasta file
    with open(fastq_file, 'r', encoding='UTF-8') as freader:
        rec_kmer = None
        unq_kmers = set()
        num_reads=0
        # Read the Strain sequence and info from the fastq file
        line_num=0
        for line in freader:
            line = line.strip()
            # process only 'n' reads from the file
            if(max_reads>0 and num_reads==max_reads):
                break

            line_num+=1
            # Ignore the quality information if the fastq file
            if(line[0] == '+' or line_num%4 == 0):
                continue

            # Start of new Strain
            if line[0] == '@':
                read_ids.append(line[1:])
                #print(f" Processing {r.split(' ')[0]} in {fasta_file}")
                # Add the t-minimizers for the earlier strain
                if(unq_kmers is not None and len(unq_kmers) > 0):
                    minimizers.append(unq_kmers)

                # Reset Heap and lists
                rec_kmer = []
                heapq.heapify(rec_kmer)
                unq_kmers = set()

            else:

                # concatenate the previous kmer sequence with the current sequence
                seq = line
                ncb_pos, ncb_cnt = check_non_canonicalbases(seq)
                #print(f' sequence length : {len(seq)}')

                i = 0
                j = 0
                last_min_pos = -1
                while i < len(seq)-k+1:

                    # Non-canonical character(s)
                    if(ncb_cnt and i <= ncb_pos[j] < (i+k)):
                        i = ncb_pos[j]  # skip kmers
                        if j < ncb_cnt-1:
                            j += 1
                    else:
                        #kmer_val = mmh3.hash(seq[i:i+k], 42, signed=False)
                        kmer_val = get64bit(mmh3.hash64(seq[i:i+k], SEED, signed=False))
                        if i < winnow_t-1 :
                            heapq.heappush(rec_kmer, (kmer_val, i))
                        else:
                            if(len(rec_kmer)>0 and kmer_val <= rec_kmer[0][0]):
                                rec_kmer = []
                                heapq.heapify(rec_kmer)

                            heapq.heappush(rec_kmer, (kmer_val, i))

                            heap_top_pos = rec_kmer[0][1]
                            if heap_top_pos != last_min_pos :
                                unq_kmers.add(rec_kmer[0][0])
                                last_min_pos = heap_top_pos

                    i += 1

                    while(len(rec_kmer)>0 and rec_kmer[0][1] <= (i-winnow_t)):
                        heapq.heappop(rec_kmer)

                num_reads+=1

        # compute the minimizers from the last read
        if(unq_kmers is not None and len(unq_kmers) > 0):
            minimizers.append(unq_kmers)
    return minimizers, read_ids

def compute_winnowing_minimizers(fasta_file: str, k: int, w: int) -> set:
    """
    Computes the minimizers based on defined winnow threshold for each Strain
    under species in the input fasta file.

    input Fasta file format:
    (odd lines) >Strain information
    (even lines) Sequence
    Returns Winnow-based minimizers for every strain under each species

    Parameters
    ----------
    fasta_file : str
    k : int
    w : int

    Returns
    -------
    set

    """
    minimizers = set([])
    # Open the fasta file
    with open(fasta_file, 'r', encoding='UTF-8') as freader:
        rec_kmer = None
        prev_seq = None
        unq_kmers = set()
        i = 0
        ic = 0
        j = 0

        # Read the Strain sequence and info from the input fasta file
        for line in freader:
            line = line.strip()

            #print('\n')

            # Start of new Strain
            if line[0] == '>':
                # Add the t-minimizers for the earlier strain
                if(unq_kmers is not None and len(unq_kmers) > 0):
                    minimizers = minimizers.union(unq_kmers)

                # Reset Heap and lists
                rec_kmer = []
                heapq.heapify(rec_kmer)
                unq_kmers = set()

                # Reset the sequence
                prev_seq = ''
                i = 0
                ic = 0

            else:

                # concatenate the previous kmer sequence with the current sequence
                seq = prev_seq+line
                ncb_pos, ncb_cnt = check_non_canonicalbases(seq)

                last_min_pos = -1
                idx = i-ic
                j = 0
                while idx < len(seq)-k+1:

                    # Non-canonical character(s)
                    if(ncb_cnt and idx <= ncb_pos[j] < (idx+k)):
                        idx = ncb_pos[j]  # skip kmers
                        if j < ncb_cnt-1:
                            j += 1
                    else:
                        #kmer_val = mmh3.hash(seq[idx:idx+k], 42, signed=False)
                        kmer_val = get64bit(mmh3.hash64(seq[idx:idx+k], SEED, signed=False))

                        if i < w-1 :
                            heapq.heappush(rec_kmer, (kmer_val, i))
                        else:
                            if(len(rec_kmer) > 0 and kmer_val <= rec_kmer[0][0]):
                                rec_kmer = []
                                heapq.heapify(rec_kmer)

                            heapq.heappush(rec_kmer, (kmer_val, i))

                            heap_top_pos = rec_kmer[0][1]
                            if heap_top_pos != last_min_pos :
                                unq_kmers.add(rec_kmer[0][0])
                                last_min_pos = heap_top_pos

                    i += 1
                    idx = i-ic

                    while(len(rec_kmer) > 0 and rec_kmer[0][1] <= (i-w)):
                        heapq.heappop(rec_kmer)

                # Store the remaining sequence for next line processing
                ic = i
                prev_seq = seq[idx:]

        # compute the minimizers from the last strain-level sequence
        if(unq_kmers is not None and len(unq_kmers) > 0):
            minimizers = minimizers.union(unq_kmers)
    return minimizers


def compute_sketch(minimizers: set, bit_size: int,
                   num_hashes: int, sketch=None) -> BloomFilter:
    """
    Create a Bloom filter based sketch using the set of minimizers

    Returns the sketch (Bloom Filter)

    Parameters
    ----------
    minimizers : set
    bit_size : int
    num_hashes : int
    sketch : BloomFilter, optional

    Returns
    -------
    BloomFilter

    """
    bit_size=bit_size*(8192) # bit_size*1024*8
    bfilter = BloomFilter(
        size=bit_size, num_hashes=num_hashes) if sketch is None else sketch
    for j in minimizers:
        bfilter.add(str(j))
    return bfilter


def write_index(sketch_nodes: list, index_file: str) -> None:
    """
    Write the index of the sketch nodes/tree into the index binary file
    The format of the index:
    (line) first 2 bytes - left node idx (even node ids),
           bit_array length - sketch/Bloom Filter size

    Parameters
    ----------
    sketch_nodes : list
    index_file : str

    Returns
    -------
    None

    """
    with open(index_file, 'wb') as ibin:
        for i in range(len(sketch_nodes)-1, -1, -1):
            lnode_idx = sketch_nodes[i].left.nidx if(
                sketch_nodes[i].left is not None) else -1
            ibin.write(lnode_idx.to_bytes(2, byteorder='little', signed=True))
            ibin.write(sketch_nodes[i].val.bit_array.tobytes())


def read_index(index_file: str, sketch_size=640, num_hashes=3) -> list:
    """
    Loads the index from the index file
    Returns the list containing left_node id and the sketch of the node
    Parameters
    ----------
    index_file : str
    sketch_size : int, optional
        DESCRIPTION. The default is 640.
    num_hashes : int, optional
        DESCRIPTION. The default is 3.

    Returns
    -------
    list

    """
    sketch_nodes = []
    sketch_size *= 1024  # in KB
    with open(index_file, 'rb') as file:
        while True:
            data = file.read(2)
            if not data:
                break
            lnode_id = int.from_bytes(data, byteorder='little', signed=True)
            bita = bitarray()
            bita.frombytes(file.read(sketch_size))
            blf = BloomFilter(sketch_size*8, num_hashes)
            blf.bit_array = bita
            sketch_nodes.append((lnode_id, blf))
    return sketch_nodes


def write_sketch_tree(sketch_nodes: list, tree_file: str) -> None:
    """
    Write the index of the sketch nodes/tree into a csv file
    The csv file has the following columns:
        1. Node index
        2. Node type (leaf nodes contain species id, intermediate nodes 'IN')
        3/4. Sketch counts (tuple)
        5. Load factor
        6. Left node index
        7. Right node index

    Parameters
    ----------
    sketch_nodes : list
    tree_file : str

    Returns
    -------
    None
    """
    # write sketch nodes from bottom-up (leaf to ROOT nodes)
    with open(tree_file, 'w', encoding='utf-8') as treef:
        for i in range(len(sketch_nodes)-1, -1, -1):
            treef.write(sketch_nodes[i].to_csv()+'\n')


def print_sketch_tree(sketch_nodes: list) -> None:
    """
    Prints the index of the sketch nodes/tree

    Parameters
    ----------
    sketch_nodes : list
    Returns
    -------
    None

    """
    # print sketch nodes from bottom-up (leaf to ROOT nodes)
    sklen = len(sketch_nodes)
    for i in range(sklen-1, -1, -1):
        print(str(sketch_nodes[i]))


def get_idx_to_id_map(tree_file: str) -> (list, dict):
    """
    Prepare Node Index to Species id, Species id to Node Index maps
    Returns Node index to Species id as list,
            Species Id to Node index as dictionary

    Parameters
    ----------
    tree_file : str

    Returns
    -------
    (list, dict)
    """
    idx_id = []
    id_idx = dict()
    with open(tree_file, 'r', encoding='utf-8') as file:
        idx=0
        for line in file:
            cols = line.split('\t')
            idx_id.append(cols[1])
            id_idx[cols[1]]=idx
            idx+=1
    return idx_id, id_idx


def get64bit(hash_tuple:tuple) -> int:
    """
    Get the 64-bit hash value by combining two 32-bit values
    Receives two 2 32-bit hash values in a tuple and
    Returns 64-bit hash value
    Parameters
    ----------
    hash_tuple : tuple
    Returns
    -------
    int
    """
    h64 = hash_tuple[1]
    h64 = (h64 << 32) + hash_tuple[0]
    return h64

def write_minimizers(minimizers:set, file:str):
    """
    Write the set of minimizers into file (for debugging purposes only)

    Parameters
    ----------
    minimizers : set
    file : str

    Returns
    -------
    None.

    """
    with open(file, 'w', encoding='UTF-8') as outf:
        for minimizer in minimizers:
            outf.write(str(minimizer)+"\n")

def query_topk_minimizers(fastq_file: str, k: int, threshold: int, max_reads=-1) -> list:
    """
    Computes the Top-k minimizers for each read in the query fastq file.
    If max_reads is set to -1, then minimizers are computed for all reads,
    else minimizers are computed upto 'max_reads' reads.

    Returns Top-k minimizers for every query read

    Parameters
    ----------
    fastq_file (query read file) : str
    k : int
    threshold : int
    max_reads : int, optional

    Returns
    -------
    list
    """
    minimizers = []
    read_ids=[]
    # Open the fasta file
    with open(fastq_file, 'r', encoding='UTF-8') as freader:
        rec_kmer = None
        unq_kmers = set()
        num_reads=0

        # Read the Strain sequence and info from the fasta file
        #for r in f.readlines():
        for line in freader:
            line = line.strip()

            # process only 'n' reads from the file
            if(max_reads>0 and num_reads==max_reads):
                break

            # Ignore the quality information if the fastq file
            if(line[0] == '+' or line[0] == '?'):
                continue

            # Start of new Strain
            if line[0] == '@':
                read_ids.append(line[1:])
                # Add the t-minimizers for the earlier strain
                if(rec_kmer is not None and len(rec_kmer) > 0):
                    minimizers.append(heapq.nsmallest(threshold, rec_kmer)
                                      if (threshold > 0) else rec_kmer)

                # Reset Heap and lists
                rec_kmer = []
                heapq.heapify(rec_kmer)
                unq_kmers = set()

            else:

                seq = line
                ncb_pos, ncb_cnt = check_non_canonicalbases(seq)

                i = 0
                j = 0
                while i < len(seq)-k+1:
                    # Non-canonical character(s)
                    if(ncb_cnt and i <= ncb_pos[j] < (i+k)):
                        i = ncb_pos[j]  # skip kmers
                        if j < ncb_cnt-1:
                            j += 1
                    else:
                        #kmer_val = mmh3.hash(seq[i:i+k], 42, signed=False)
                        kmer_val = get64bit(mmh3.hash64(seq[i:i+k], SEED, signed=False))
                        if kmer_val not in unq_kmers:
                            unq_kmers.add(kmer_val)
                            heapq.heappush(rec_kmer, kmer_val)
                    i += 1

                num_reads+=1

        # compute the minimizers from the last strain-level sequence
        if len(rec_kmer) > 0:
            minimizers.append(heapq.nsmallest(threshold, rec_kmer)
                              if (threshold > 0) else rec_kmer)
    return minimizers, read_ids

def compute_topk_minimizers(fasta_file: str, k: int, threshold: int) -> set:
    """


    Parameters
    ----------
    fasta_file : str
        DESCRIPTION.
    k : int
        DESCRIPTION.
    threshold : int
        DESCRIPTION.

    Returns
    -------
    set
        DESCRIPTION.

    """
    minimizers = set([])
    # Open the fasta file
    with open(fasta_file, 'r', encoding='UTF-8') as freader:
        rec_kmer = None
        prev_seq = None
        unq_kmers = set()

        # Read the Strain sequence and info from the fasta file
        for line in freader:
            line = line.strip()

            # Start of new Strain
            if line[0] == '>':
                #print(f" Processing {r.split(' ')[0]} in {fasta_file}")
                # Add the t-minimizers for the earlier strain
                if(rec_kmer is not None and len(rec_kmer) > 0):
                    minimizers = minimizers.union(heapq.nsmallest(
                        threshold, rec_kmer) if (threshold > 0) else rec_kmer)

                # Reset Heap and lists
                rec_kmer = []
                heapq.heapify(rec_kmer)
                unq_kmers = set()

                # Reset the sequence
                prev_seq = ''

            else:

                # concatenate the previous kmer sequence with the current sequence
                seq = prev_seq+line
                ncb_pos, ncb_cnt = check_non_canonicalbases(seq)

                i = 0
                j = 0
                while i < len(seq)-k+1:
                    # Non-canonical character(s)
                    if(ncb_cnt and i <= ncb_pos[j] < (i+k)):
                        i = ncb_pos[j]  # skip kmers
                        if j < ncb_cnt-1:
                            j += 1
                    else:
                        #kmer_val = mmh3.hash(seq[i:i+k], 42, signed=False)
                        kmer_val = get64bit(mmh3.hash64(seq[i:i+k], SEED, signed=False))
                        if kmer_val not in unq_kmers:
                            unq_kmers.add(kmer_val)
                            heapq.heappush(rec_kmer, kmer_val)
                    i += 1

                # Store the remaining sequence for next line processing
                prev_seq = seq[i:]

        # compute the minimizers from the last strain-level sequence
        if len(rec_kmer) > 0:
            minimizers = minimizers.union(heapq.nsmallest(
                threshold, rec_kmer) if (threshold > 0) else rec_kmer)

    return minimizers
