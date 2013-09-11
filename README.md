Glust
===========
Gclust (Genome sequence clustering), a program for clustering the rapid growth of complete or draft genome sequences. Using a sparse suffix array algorithm and a genomic-distance identity taking into account the criteria of diversity, which is based on extension DNA maximal exact matches (MEM), Gclust creates clusters under the given set of genome sequences and extension MEM identity. It takes less than 7 hours for the clustering of the 1560 complete microbial genome sequences with average 3.4MB length on Intel(R) Xeon(R) CPU 2.27GHz with 8 threads parallel computing. It offers the possibility of clustering the rapid growth of complete or draft microbial genomes in the future.

Gclust is special designed for genome sized sequences clustering and introduced one kind of genomic-distance identity taking into account the criteria of diversity. The fast sparse suffix array construction algorithm was used in finding MEMs between query genome sequence and representative genome sequences. The dynamic programming extension of MEMs is also supported for genome sequence identity computing. Our implementation supports multithreads parallel computing.

Gclust was written in C++ and uses the SeqAn library and the libdivsufsort library. It is currently maintained by Dr. Beifang Niu (bniu@genome.wustl.edu) in the Dr. Ding's group (http://genome.wustl.edu/people/groups/detail/ding-lab).

Usage
-----

        Version 0.5
        Usage:  gclust [options] <genome-file>

Options:

Usage: ./gclust [options] <genomes-file>

Options:

    -minlen        set the minimum length of a exact match, the default value is 20
    -both          compute forward and reverse complement matches, the default value is only forward strand
    -nuc           match only the characters a, c, g, or t, the defult value is yes
    -sparse        set the step of sparse suffix array, the default value is 1
    -threads       set the number of threads to use, the default value is 1
    -block         set the block size for one time clustering, the default value is 100 (MB)
    -nblock        set the blocks number loaded one time for remaining genomes alignment, the default value is 2
    -loadall       loading total genomes one time, the default value is yes
    -rebuild       rebuild suffix array after clustering into one block, the default value is yes

Clustering cutoff:

    -memiden       MEMs identity for clustering, the default value is 90 ( 90% identity )

MEM extension options:

    -ext       Extension options, 0: No extension, 1: Gapped extension, 2: Ungapped extension, the default is 1
    -mas       Reward for a nucleotide match, the default value is 1
    -umas      Penalty for a nucleotide mismatch, the default value is -1
    -gapo      Cost to open a gap, the default value is -1
    -gape      Cost to extend a gap, the default value is -1
    -drops     X dropoff value for extension, the default value is 1

Example usage:

    ./gclust -minlen 20 -both -nuc -threads 8 -ext 1 genome.fa > clusering.out

Find all gapped extension MEMs on forward and reverse strands of length 20 or greater, matching only a, c, t, or g using 8 threads parallel computing.


Install
-------
    make
    make install

xxx

