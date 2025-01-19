The Department of Mathematics,  Amrita Vishwa Vidyapeetham, Coimbatore, organized “ACM School on Data Structures and Algorithms for Strings, with Applications to Search Engines and Computational Biology”. This assignment is a part of the coursework of Winter School.


## Querying the Microbial sequences with compact and accurate index

### Download and Install
*Please ensure atleast 10 GB of hard disk space.*

1. Download the assignment.tar.gz and prepare.sh (or) prepare.bat
2. Run installation using prepare.sh
   ```Shell
   chmod 755 prepare.sh
   ./prepare.sh
   ```
3. The project folder structure should be as follows:

- source code files
- data (data directory)
  - input
  - query
  - index
  - results

### 1. Code
The Code repository contains the following python source files:
1. index.py : Reads the input sequences (.fasta) and creates the microbial sequence index for all species
2. query.py : Reads the input query sequences (.fastq), queries the index, computes the Top-k metrics across the species.
3. helper.py : Utility functions that enable construction of sketches (including winnow minimizers with bloom filters), read and write of index.
4. Node.py : Node class for tree representation of index
5. BloomFilter.py : BloomFilter class with add, find, intersection and union operations.

### 2. Index
The folder contains several indexes created using different parameters. Each index (tree) is represented and stored in binary and csv format.

### 3. Query
The folder contains query read sequences in the fastq format for 10 species.
Format : Fastq file denoting species id. Lines starting with @ contains read information
Lines starting with + and every fourth line contains read quality information
Every second line contains the actual read
```
Line 1: @readsim-datasets/1771-readsim_000000001_L000010000:000000002-000010001:F
Line 2: GACAGCGGACCCCGACCCACCTTTCGTGCCGTCTGGAACACTGTCGTAGCCGAACTCAACGCCG...
Line 3: +
Line 4: MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM2MMM .. 
```

- 10 species
- Reads per species : 1000
- Seq. lengths per read : 10000
- Total disk space (uncompressed) : 200.6 MB
- Sequences are only forward strand
    - Contains Non-canonical bases (BKNDYWMVRH)
    - Contains Mutations, Sequencing noise/artefacts.

### 4. Input
The folder contains input read sequences in the fasta format for 1024 species.

Format : Fasta file denoting species id. Lines starting with > contains Strain information
Every other line contains raw sequences till the next strain information line.
```
Line 1: >NZ_CP040018.1 Arthrobacter sp. 24S4-2 chromosome, complete genome
Line 2: ATGACAGTAGACGAAGCCAACCACGCGAATACTGTCGGAAGTTCCTGGCGGCGGGTTGTGAGCCTCCTGGAGCAGGACCA
Line 3: CCGGGTTTCACCCCGGCAGCGCGGCTTCGTAATCCTCGCCCAGGCGCAGGGACTGATCGGATCCACCCTCTTGGTGGCCG
```

- 1024 species, 1914 strains
- Strains per species : 1-25
- Seq. lengths across strains : 5.32 – 15.8 mn
- Total disk space (uncompressed) : 8.6 GB
- Sequences are only forward strand
    - Contains Non-canonical bases (BKNDYWMVRH)
    - Contains Mutations, Sequencing noise/artefacts.

## Indexing phase
The Microbial index creation involves:
1. Reading the input sequence fasta files across all species
2. Computing sketches for each of the species
3. Building a index (tree) for the all the sketches
4. Writing the index

Indexes can created by executing the code below. 
- The default parameters :
(Minimizer based sketches (-m) , 
kmer size (-k) : 21, 
winnow length (-w) : 50, 
sketch size (-s) : 640 KB, 
number of hash functions (-n) : 3)

```python
    python3 -u index.py -d <data_dir>
```

## Querying and metrics computation phase
Querying microbial index involves:
1. Reading the specified index from the disk
2. Computing the sketches for each of the query species.
3. Querying the index with the computed query sketches
4. Calculate and write the Top-K metrics.

Querying the index can be performed by executing the code below.
- The default parameters :
(Minimizer based sketches (-m) , 
similarity threshold (-s) : 0.7, 
maximum reads/sequences (-n) : -1 (all reads))

```python
    python3 -u query.py -d <data_dir> -i <index_file>
```
