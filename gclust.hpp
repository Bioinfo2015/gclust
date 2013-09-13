/*
 * gclust.hpp for gclust
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <cctype>
#include "fasta.hpp"
#include "paraSA.hpp"

using namespace std;

void usage(string prog);

int K = 1; // Note: Using sparse suffix array for larger chunk size.
int Nchunk = 2; // load one time for remaining genomes clustering.
int chunk = 100; // block size for clustering part by part.
int min_len = 20; // Default minimum exact match length.
int MEMiden = 90; // Default identity cutoff.
int total_threads = 1; // Threads number.

// MEM extension parameters
int ext = 1; // noextension, gap or ungap extension
int mas = 1; // Match score
int umas = -1; // Mismatch cost
int gapo = -1; // Gap open penalty
int gape = -1; // Gap extension penalty
int drops = 1; // Maximum score drop

bool rev_comp = false;
bool nucleotides_only = true;
bool rebuild = false; // Rebuild suffix array into one part.
bool loadall = false; // load all genomes one time, need more memory.

paraSA *sa, *saa; // Suffix array.

vector< Genome > refseqs, allrefseqs; // Part genomes and total part genomes.
vector< GenomeClustInfo > totalgenomes; // Total genomes.
vector< vector<match_t> > matchlist; // Parallel buffer for match_t.
vector< vector<mumi_unit> > mumilist; // parallel buffer for mumi_unit.

//Multithreads parallel parameters passing
struct threads_arg { 
    int skip;
    int skip0; // Zebra-style distribution
    long chunk;
    long begin; // Begin alignment position
    bool part; // If internal part clustering
};

// Note: just test distances between genomes.
void testDistanceBgenomes(vector<GenomeClustInfo> &totalgenomes) {
    for (long i=0; i<(long)totalgenomes.size(); i++) {
        if (totalgenomes[i].clusters.size()>0) {
            for (long j=0; j<(long)totalgenomes[i].clusters.size(); j++) {
                hit thit;
                thit = totalgenomes[i].clusters[j];
                cerr<<thit.id<<"\t"<<thit.identity<<endl;
            }
        }
    }
}

// Note: output clustering information as cd-hit format.
void outputClusteringInfoSimple(vector<GenomeClustInfo> &totalgenomes) {
    long clusters = 0;
    for (long i=0; i<(long)totalgenomes.size(); i++) {
        if (totalgenomes[i].rep) {
            long clusterunit=0;
            // Output representive.
            cout << ">Cluster " << clusters << endl;
            cout << clusterunit
                 <<"\t" << totalgenomes[i].size
                 << "nt, >" << totalgenomes[i].descript 
                 << "... *" << endl;
            if (totalgenomes[i].clusterunits.size()>0) {
                for (long j=0;j<(long)totalgenomes[i].clusterunits.size();j++) {
                    clusterunit++;
                    cout << clusterunit << "\t" 
                         << totalgenomes[totalgenomes[i].clusterunits[j].id].size<<"nt, >"
                         << totalgenomes[totalgenomes[i].clusterunits[j].id].descript 
                         << "... at " << totalgenomes[i].clusterunits[j].strand 
                         << "/" << totalgenomes[i].clusterunits[j].identity
                         << endl;
                }
            }
            clusters++;
        }
    }
    cerr << "Total clusters: " << clusters << endl;
}

// Note: collect clustering information.
void getClusteringInfoOnepart(vector<GenomeClustInfo> &totalgenomes, long begin, long chunk, bool inpart, bool &clusterhit) {
    long b,e;
    clusterhit = false;
    b = begin;
    e = begin+chunk;
    for (long i=b; i<e; i++) {
        if (!totalgenomes[i].rep){ continue; }
        if (totalgenomes[i].clusters.size()>0) {
            for (long j=0; j<(long)totalgenomes[i].clusters.size(); j++) {
                hit thit;
                thit = totalgenomes[i].clusters[j];
                if (totalgenomes[thit.id].rep) {
                    hit thit0;
                    thit0.id = i;
                    thit0.identity = thit.identity;
                    thit0.strand = thit.strand;
                    totalgenomes[thit.id].clusterunits.push_back(thit0);
                    totalgenomes[i].rep = false;
                    clusterhit = true;
                    break;
                }
            }
            // Note: important clear
            totalgenomes[i].clusters.clear();
        }
    }
}

// Note: one genome as reference (internal part).
void *single_thread(void *arg_) {
    Genome tg;
    threads_arg *arg = (threads_arg *)arg_;
    // Match information container.
    vector<match_t> &matches = matchlist[arg->skip0];
    // Mem index container.
    vector<mumi_unit> &mumis = mumilist[arg->skip0];
    
    long seq_cnt = 0;
    long beginclust = arg->begin;
    long chunk = arg->chunk;
    long sizeadd = 0;
    long edge = long(refseqs.size()-1);
    bool ifhit = false;
    bool ispart = arg->part;
    double cutoff=(double)MEMiden/100;
    string *P = new string; 
    edge = long(refseqs.size()-1);
    
    while (1) {
        if ( seq_cnt > edge ) { break; }
        if ( arg->skip0 == 0 ) {
            if (seq_cnt % 100 ==0) { cerr <<"...... "<<seq_cnt<<" done"<< endl; }
        }
        // paralle part
        if ( seq_cnt % arg->skip == arg->skip0 ) {
            ifhit = false;
            tg = refseqs[seq_cnt];
            if ( totalgenomes[tg.id].rep ) {
                *P = tg.cont;
                // Filter 'n'
                if (nucleotides_only) { filter_n(*P); }
                // 100% ?
                if (MEMiden == 100) {
                    saa->MEMperfect(*P, matches, tg.size, tg.id);
                } else { saa->MEM(*P, matches, min_len, tg.id); }
                sizeadd += saa->load_match_info(tg.id, matches, mumis, true, tg.size);
                matches.clear();
                if ((double)sizeadd/tg.size >= cutoff) {
                    ifhit = ComputeMemIdentity(totalgenomes, allrefseqs, mumis, beginclust, tg.id, MEMiden, ispart, chunk, '+', ext,  mas, umas, gapo, gape, drops);
                }
                mumis.clear();
                sizeadd = 0;
                if ((ispart) || (!ifhit)) {
                    if(rev_comp) {
                        reverse_complement(*P, nucleotides_only);
                        // 100% ?
                        if (MEMiden == 100) {
                            saa->MEMperfect(*P, matches, tg.size, tg.id);
                        } else {  saa->MEM(*P, matches, min_len, tg.id); }
                        // Loading match information - strand.
                        sizeadd += saa->load_match_info(tg.id, matches, mumis, true, tg.size);
                        matches.clear();
                        if ((double)sizeadd/tg.size >= cutoff) {
                            ComputeMemIdentity(totalgenomes, allrefseqs, mumis, beginclust, tg.id, MEMiden, ispart, chunk, '-', ext, mas, umas, gapo, gape, drops);
                        }
                        sizeadd = 0;
                        mumis.clear();
                    }
                }
            }
        }
        seq_cnt++;
        delete P; 
        P = new string;
    }

    delete P;
    pthread_exit(NULL);

}


void usage( string prog )  {
    cerr << "Gclust is a clustering program for genome, draft assembly contigs,";
    cerr << " which algorithm is based on all Maximal Exact Matches(MEMs) between genome sequences." << endl;
    cerr << endl;
    cerr << "Usage: " << prog << " [options] <genomes-file> " << endl;
    cerr << endl;
    cerr << "Options:" << endl;
    cerr << endl;
    cerr << "-minlen        set the minimum length of a exact match, the default value is 20" << endl;
    cerr << "-both          compute forward and reverse complement matches, the default value is only forward strand" << endl;
    cerr << "-nuc           match only the characters a, c, g, or t, the defult value is yes" << endl;
    cerr << "-sparse        set the step of sparse suffix array, the default value is 1" << endl;
    cerr << "-threads       set the number of threads to use, the default value is 1" << endl;
    cerr << "-block         set the block size for one time clustering, the default value is 100 (MB)" << endl;
    cerr << "-nblock        set the blocks number loaded one time for remaining genomes alignment, the default value is 2" << endl;
    cerr << "-loadall       loading total genomes one time, the default value is yes" << endl;
    cerr << "-rebuild       rebuild suffix array after clustering into one block, the default value is yes" << endl;
    cerr << endl;
    cerr << "Clustering cutoff:" << endl;
    cerr << endl;
    cerr << "-memiden       MEMs identity for clustering, the default value is 90 ( 90% identity )" << endl;
    cerr << endl;
    cerr << "MEM extension options:" << endl;
    cerr << endl;
    cerr << "-ext       Extension options, 0: No extension, 1: Gapped extension, 2: Ungapped extension, the default is 1 " << endl;
    cerr << "-mas       Reward for a nucleotide match, the default value is 1 " << endl;
    cerr << "-umas      Penalty for a nucleotide mismatch, the default value is -1 " << endl;
    cerr << "-gapo      Cost to open a gap, the default value is -1 " << endl;
    cerr << "-gape      Cost to extend a gap, the default value is -1 " << endl;
    cerr << "-drops     X dropoff value for extension, the default value is 1 " << endl;
    cerr << endl;
    cerr << "Example usage:" << endl;
    cerr << endl;
    cerr << "./gclust -minlen 20 -both -nuc -threads 8 -ext 1 genome.fa > clusering.out" << endl;
    cerr << endl;
    cerr << "Find all gapped extension MEMs on forward and reverse strands" << endl;
    cerr << "of length 20 or greater, matching only a, c, t, or g" << endl;
    cerr << "using 8 threads parallel computing." << endl;
    cerr << endl;
    exit(1);
}

