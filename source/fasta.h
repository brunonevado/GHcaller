//
// fasta.h
//  NGHcaller
//
//  Copyright (c) 2013 Bruno Nevado, GNU license.

#ifndef __ngh_fasta__
#define __ngh_fasta__

#include <iostream>
#include <vector>
#include <string>

#include "common.h"


struct fasta_bed {
    std::vector <unsigned int > rstarts, rends;
    std::vector <std::string> rnames;
    std::string infile, chr;
};


class fasta {
    std::vector <std::string> matrix;
    std::vector <std::string> names;
    std::string infile;
    fasta_bed regions;
    void write_to_file_haplos ( std::string out, int start0, int end0 );
    void write_to_file_iupac ( std::string out, bool rooted , int start0, int end0 );
public:
    fasta( int num_inds, int len = 0 );
    unsigned int num_lines () const {return int ( matrix.size() );}
    unsigned int num_bases () const {return int ( matrix[0].size() );}
    void set_num_inds ( unsigned int value ) { matrix.resize(value); }
    void set_all_names ( std::vector<std::string> in );
    void read_fasta_file ( std::string in ) ;
    std::string name_at ( int ind0  ) { return names.at(ind0); };
    void write ( std::string out, int wlen, bool rooted, int iupac, bool verb  );
    void write_bed ( std::string out, int wlen, bool rooted, int iupac, bool concat, bool verb  );
    void build_from_ngh( unsigned int site, std::string col_to_add );
    void append_seq (fasta in );
    void check_bed (std::string file_name);
};

#endif