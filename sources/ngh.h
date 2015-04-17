//
//  ngh.h
//  GHcaller
//
//  Copyright (c) 2013 Bruno Nevado, GNU license.

#ifndef __ngh__
#define __ngh__

#include <iostream>
#include <vector>
#include <string>

#include "common.h"

struct ngh_parse_return {
    std::string chr;
    unsigned int site;
    char ref;
    unsigned int ninds;
    std::vector <std::string> reads;
    std::vector <std::string> quals;
    std::vector <int> raw_deps;
};

struct res_snp_call {
    std::string called_genotypes;
    std::vector < std::vector <float>> gls; // AA, AC, AG, AT, CC, CG, CT, GG, GT, TT, for each ind
};



class ngh {
private:
    std::vector <unsigned int> mindeps, maxdeps, platforms, baseqs;
    unsigned int num_inds, minreads;
    std::vector <char> DNAbases;
    std::vector < std::vector < char > > possible_genotypes;
    int winlen, iupac;  // options for fasta output
    int lrt_genotypes ( std::vector <float> gls );
    int lrt_haplotypes ( std::vector <float> gls );
    float p_err, chivalue;
    std::vector <float> calc_genos_probs ( std::string reads );
    float calc_homo_p ( int n1, int n2, int n3, int n4 );
    float calc_hetero_p ( int n1, int n2, int n3, int n4 );
    bool do_haplos;
public:
    ngh( std::string mindep, std::string maxdep, std::string platform, std::string baseq,
        int minreads, float err, float chi, int wlen, int iupac, bool haplos );
    ngh_parse_return parse_mpileup( std::string  line );
    res_snp_call roesti_snpcall ( ngh_parse_return in );
    std::vector < std::string> info_to_vector();
    void set_num_inds (  unsigned int );
};




#endif
