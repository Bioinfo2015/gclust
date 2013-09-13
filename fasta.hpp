
/*
 * fasta.hpp for gclust
 * Copyright (c) 2010-2011 Beifang Niu All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */


#ifndef __FASTA_HPP__
#define __FASTA_HPP__

#include <string>
#include <vector>

using namespace std;

const long MAX_GENOME = 1000000000;
const long MAX_THREADCONTAINER = 10000000;
const long MAX_PARTNUMBER = 500;
const long MAX_PARTNUMBERFORPERFECT = 40000;
const long PART_BASE = 1000000; // 1MB

// Clustering information
struct hit {
    long id;
    double identity;
    char strand;
};

// For single genome search
struct Genome {
    long size; //genome length
    long id; //index
    string descript; // genome name
    string cont; // genome content
};

// Total genomes info container
struct GenomeClustInfo {
    long size; //genome length
    long id; //index
    bool rep; // if it is representive of clustering
    double mumindex; //another index
    string descript; //genome name
    vector<hit> clusters, clusterunits;
};

void filter_n( string &seq_rc );
// Reverse complement sequence
void reverse_complement( string &seq_rc, bool nucleotides_only);
// Trim sequence
void trim( string &line, long &start, long &end );
void load_fasta( string filename, string &S, vector<string> &descr, vector<long> &startpos );
void load_part_genomes( string filename, vector<Genome> &partgenomes, vector<GenomeClustInfo> &totalgenomes, long previous, long number );
void load_part_genomes_internal( string filename, vector<Genome> &partgenomes, vector<GenomeClustInfo> &totalgenomes, long previous, int &number, long totalsize, bool &ifend, int memiden );
void load_part_genomes_all( string filename, vector<Genome> &partgenomes );
// Load part genomes into memory between parts
void load_part_genomes_mem( vector<Genome> &allpartgenomes, vector<Genome> &partgenomes, vector<GenomeClustInfo> &totalgenomes, long previous, long number );
// Load part genomes into memory
void load_part_genomes_internal_mem( vector<Genome> &allpartgenomes, vector<Genome> &partgenomes, vector<GenomeClustInfo> &totalgenomes, long previous, int &number, long totalsize, bool &ifend, int memiden );
void test_part( vector<Genome> &partgenomes );
// Make block for Suffix Array construction
void make_block_ref( vector<Genome> &partgenomes, string &S, vector<GenomeClustInfo> &totalgenomes, vector<long> &descr, vector<long> &startpos );
// Load total genomes
void load_total_genomes( string filename, vector<GenomeClustInfo> &totalgenomes );

#endif // __FASTA_HPP__

