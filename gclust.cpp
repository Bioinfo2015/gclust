
/*
 * gclust.cpp for gclust
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
#include "gclust.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    // Version notice
    cerr <<"\nGclust version 0.5\n"<< endl;
    // Collect parameters from the command line.
    while (1) {
        static struct option long_options[] =
        { 
            {"minlen", 1, 0, 0}, // 0
            {"both", 0, 0, 0}, // 1
            {"nuc", 0, 0, 0}, // 2
            {"threads", 1, 0, 0}, // 3
            {"block", 1, 0, 0}, //4
            {"memiden", 1, 0, 0}, //5
            {"nblock", 1, 0, 0}, //6
            {"loadall", 0, 0, 0}, //7
            {"rebuild", 0, 0, 0}, //8
            // Seed extension part
            {"mas", 1, 0, 0}, //9
            {"umas", 1, 0, 0}, //10
            {"gapo", 1, 0, 0}, //11
            {"gape", 1, 0, 0}, //12
            {"drops", 1, 0, 0}, //13
            {"ext", 1, 0, 0}, //14
            // Sparse step of suffix array
            {"sparse", 1, 0, 0,}, //15
            {0, 0, 0, 0}
        };
        int longindex = -1;
        int c = getopt_long_only(argc, argv, "", long_options, &longindex);
        if (c == -1) break; 
        else if (c == '?') { 
            cerr << "Invalid parameters." << endl;
            usage(argv[0]);
        } else {
            // Branch on long options
            switch(longindex) { 
                case 0: min_len = atol(optarg); break;
                case 1: rev_comp = true;	break;
                case 2: nucleotides_only = true; break;
                case 3: total_threads = atoi(optarg) ; break;
                case 4: chunk = atoi(optarg) ; break;
                case 5: MEMiden = atoi(optarg) ; break;
                case 6: Nchunk = atoi(optarg) ; break;
                case 7: loadall = true ; break;
                case 8: rebuild = true ; break;
                // Seed extension part
                case 9: mas = atoi(optarg) ; break;
                case 10: umas = atoi(optarg) ; break;
                case 11: gapo = atoi(optarg) ; break;
                case 12: gape = atoi(optarg) ; break;
                case 13: drops = atoi(optarg) ; break;
                case 14: ext = atoi(optarg) ; break;
                // Sparse step of suffix array
                case 15: K = atoi(optarg) ; break;
                default: break; 
            }
        }
    }
    // Only using all maximal matches for clustering
    if (argc - optind != 1) usage(argv[0]);
    if (total_threads <= 0) { 
        cerr << "invalid number of threads specified" << endl; 
        exit(1); 
    }
    // no extension when 100% match
    if (MEMiden == 100) { ext = 0; }
    // Allocate memory for multithreads
    for (int i=0; i<total_threads; i++) {
        vector<match_t> matches;
        vector<mumi_unit> mumis;
        matchlist.push_back(matches);
        mumilist.push_back(mumis);
        matchlist[i].reserve(MAX_THREADCONTAINER);
        mumilist[i].reserve(MAX_THREADCONTAINER);
    }
    // Genome file
    string ref_fasta = argv[optind]; 
    // Load total genomes part information
    load_total_genomes(ref_fasta, totalgenomes);
    // Load total part genomes one time
    if (loadall) load_part_genomes_all(ref_fasta, allrefseqs);
    
    vector<long> refdescr;
    vector<long> startpos;
    long chunksize;
    long dchunk;
    long genomes;
    long begin = 0;
    long cbegin = 0;
    bool ifend=false;
    bool clusterhit=false;
    Genome tg; 
    string ref;
    // set chunk size for clustering chunk by chunk
    chunksize=(long)chunk*PART_BASE;
    // Note: Take a fixed value of total genomes number
    // Main clustering loop
    while (1) {
        if (loadall) { 
            load_part_genomes_internal_mem(allrefseqs, refseqs, totalgenomes, cbegin, chunk, chunksize, ifend, MEMiden);
        } else {
            load_part_genomes_internal(ref_fasta, refseqs, totalgenomes, cbegin, chunk, chunksize, ifend, MEMiden);
        }
        cerr <<"\nLoad genomes: "<<refseqs.size()<< endl;
        cerr <<"\nchunk: "<<chunk<< endl;
        refdescr.clear(); startpos.clear();
        // Clear container
        ref = ""; 
        // Make part suffix array
        make_block_ref(refseqs, ref, totalgenomes, refdescr, startpos);
        cerr <<"Creating suffix array ......\n"<< endl;
        saa = new paraSA(ref, refdescr, startpos, true, K);
        cerr <<"\nFinished creating suffix array ......\n"<< endl;
        genomes = refseqs.size();
        // Part internal clustering || parallel part
        pthread_attr_t attr;  pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        vector<threads_arg> args(total_threads);
        vector<pthread_t> thread_ids(total_threads);
        // Initialize additional thread data
        for (int i=0; i<total_threads; i++) {
            args[i].skip = total_threads;
            args[i].skip0 = i;
            args[i].begin = cbegin;
            args[i].part = true;
            args[i].chunk = chunk; 
        }
        // Create joinable threads to find MEMs
        for (int i=0; i<total_threads; i++) pthread_create(&thread_ids[i], &attr, single_thread, (void *)&args[i]);
        // Wait for all threads to terminate
        for (int i=0; i<total_threads; i++) pthread_join(thread_ids[i], NULL);
        // Collect clustering information into one chunk
        getClusteringInfoOnepart(totalgenomes, cbegin, chunk, true, clusterhit);
        if (ifend) { break;} // Finished
        if ( (rebuild)&&(clusterhit) ) {
            delete saa; ref="";
            refdescr.clear(); 
            startpos.clear();
            // Make part suffix array
            make_block_ref(refseqs, ref, totalgenomes, refdescr, startpos);
            cerr <<"Creating suffix array ......\n"<< endl;
            saa = new paraSA(ref, refdescr, startpos, true, K);
            cerr <<"\nFinished creating suffix array ......\n"<< endl;
        }
        
        refseqs.clear();
        begin = cbegin+chunk;
        if (loadall) { dchunk=(long)allrefseqs.size(); }else{ dchunk=(long)(chunk*Nchunk); }
        // Make alignment for last genomes
        while (1) {
            cerr <<"\n=================="<< endl;
            cerr <<"begin alignment "<<begin<<"\n"<< endl;
            if (loadall) { 
                load_part_genomes_mem(allrefseqs, refseqs, totalgenomes, begin, dchunk);
            } else {
                load_part_genomes(ref_fasta, refseqs, totalgenomes, begin, dchunk);
            }
            // Parallel part
            pthread_attr_t attr;  pthread_attr_init(&attr);
            pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            vector<threads_arg> args(total_threads);
            vector<pthread_t> thread_ids(total_threads);
            // Initialize additional thread data
            for (int i = 0; i < total_threads; i++) {
                args[i].skip = total_threads;
                args[i].skip0 = i;
                args[i].begin = cbegin;
                args[i].part = false;
                args[i].chunk = chunk;
            }
            // Create joinable threads to find MEMs
            for (int i = 0; i < total_threads; i++) pthread_create(&thread_ids[i], &attr, single_thread, (void *)&args[i]);
            // Wait for all threads to terminate
            for (int i = 0; i < total_threads; i++) pthread_join(thread_ids[i], NULL);
            if ((long)refseqs.size()<dchunk){ 
                getClusteringInfoOnepart(totalgenomes, begin, (long)refseqs.size(), false, clusterhit); 
            } else {
                getClusteringInfoOnepart(totalgenomes, begin, dchunk, false, clusterhit);
            }
            if ((long)refseqs.size()<dchunk) { break; }
            begin += dchunk;
            refseqs.clear();
        }
        delete saa;
        refseqs.clear();
        cbegin = cbegin + chunk;
    }//end while(1)
    //testDistanceBgenomes(totalgenomes); 
    // Collect clustering information
    cerr <<"\n==========================="<< endl; 
    cerr <<"Output clustering information ......\n"<< endl;
    // Output with CD-HIT format
    outputClusteringInfoSimple(totalgenomes);
}

