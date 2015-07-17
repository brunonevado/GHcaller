//
//  fasta.cpp
//  GHcaller
//
//  Copyright (c) 2013 Bruno Nevado, GNU license.

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <locale>

#include "fasta.h"

// FASTA OBJECT CONSTRUCTOR

fasta::fasta ( int num_inds, int len ){
    if(len == 0){
        matrix.reserve(num_inds);
    }
    else{
        matrix.resize(num_inds);
        
        for(int i = 0; i < num_inds; i++)
            matrix.at(i).resize(len);
    }
}


// FASTA READ FASTA FILE

void fasta::read_fasta_file( std::string infas){
    int cind = 0;
    matrix.clear();
    std::locale loc;
    std::string line;
    std::ifstream infile_fas (infas.c_str());
    if (infile_fas.is_open())
    {
        infile = infas;
        while ( ! infile_fas.eof() )
        {
            getline(infile_fas, line);
            
            if( line[0] == '>' ){
                cind++;
                matrix.resize(cind);
                std::string name = line.substr(1);
                names.push_back(name);
            }
            else {
                matrix.at(cind-1).append( line );
            }
        }
    }
    
    else {
        std::cerr << "ERROR (read_fasta_file): Unable to open infile " << infas << "\n" ;
        exit(1);
    }
    infile_fas.close();
    
    for (unsigned int l = 0; l < matrix.size(); l++) {
        for (unsigned s = 0; s < matrix.at(l).length(); s++) {
            matrix.at(l).at(s) = tolower(matrix.at(l).at(s), loc);
        }
    }
    
}


// WRITE ALIGNMENT TO FILE

void fasta::write(std::string out, int wlen, bool rooted, int iupac, bool verb ){
    
    int nwindows = (wlen != 0) ? this->num_bases()/wlen : 1;
    if (wlen != 0 && this->num_bases() % wlen != 0) {
        nwindows++;
    }
    
    for (int n = 0; n < nwindows; n++) {
        int start = n*wlen+1;
        unsigned int end = (n*wlen+ wlen <= this->num_bases() && wlen != 0 )? n*wlen+ wlen : this->num_bases();
        std::stringstream outfile;
        outfile << out;
        if( nwindows > 1){
            outfile.str("");
            outfile << out << "." << start << "_" << end << ".fas";
        }
        if(verb)
            std::clog << "<GHcaller> Writing window " << n+1 << " to file " << outfile.str() <<", region:" << start << "-" << end << std::endl;
        
        
        if (iupac == 1) {
            this->write_to_file_iupac(outfile.str(), rooted, start - 1, end -1 );
        }
        else if (iupac == 0){
            this->write_to_file_haplos(outfile.str(), start - 1, end -1 );
        }
        
    }
    
}

// WRITE BED
void fasta::write_bed(std::string out, int wlen, bool rooted, int iupac, bool concat, bool verb ){
 
    if(!concat){
        for (unsigned int n = 0; n < regions.rstarts.size(); n++) {
            
            std::stringstream outfile;
            outfile << out;
            if( regions.rstarts.size() > 1){
                outfile.str("");
                outfile << out << "." << regions.rnames.at(n) << ".fas";
            }
            if(verb)
                std::clog << "<GHcaller> Writing region " << regions.rnames.at(n) << " to file " << outfile.str() <<", region:" << regions.rstarts.at(n) << "-" << regions.rends.at(n) << std::endl;
            
            
            if (iupac == 1) {
                this->write_to_file_iupac(outfile.str(), rooted, regions.rstarts.at(n) -1 , regions.rends.at(n) -1 );
            }
            else if (iupac == 0){
                this->write_to_file_haplos(outfile.str(), regions.rstarts.at(n) -1 , regions.rends.at(n) -1  );
            }
        }
    }
    else{
        fasta aBedFasta(1);
        aBedFasta.set_num_inds( this->num_lines() );
        
        for (unsigned int iname = 0; iname < names.size(); iname++) {
            aBedFasta.names.push_back(names.at(iname));
        }
        
        
        for (unsigned int n = 0; n < regions.rstarts.size(); n++) {
            for (unsigned int i = 0; i < matrix.size(); i++) {
                aBedFasta.matrix.at(i).append( matrix.at(i).substr( regions.rstarts.at(n)-1,  regions.rends.at(n) - regions.rstarts.at(n) + 1 ) );
            }
            
        }
        
        aBedFasta.write(out, wlen, rooted, iupac, false);
         
        if(verb){
            std::clog << "<GHcaller> Fasta file written to " << out << " with all regions defined in bed file " << regions.infile << std::endl;
        }
        
    }
    
    
}

void fasta::write_to_file_haplos( std::string out, int start0, int end0 ){
    std::locale loc;
    std::ofstream outputFile;
    // set haplo names
    for (unsigned int n = 1; n < names.size(); n += 2) {
        names.at(n-1).append("_1");
        names.at(n).append("_2");
    }
    
    outputFile.open(out.c_str());
    
    if( !outputFile.is_open() ){
        std::cerr << "ERROR (write_to_file): unable to open for output file " << out << std::endl;
        exit(1);
    }
    
   
    for( unsigned int i = 0; i < matrix.size(); i++){
        outputFile << ">" << names[i] << "\n";
        
        for ( int site = start0; site <= end0; site++) {
            outputFile << toupper(matrix[i][site],loc);
        }
        outputFile << std::endl;
    }
    
    outputFile.close();
    
}

void fasta::write_to_file_iupac ( std::string out, bool rooted, int start0, int end0 ){
    // if rooted, wont make IUAPC of last line
    std::ofstream outputFile;
    std::locale loc;
    
    outputFile.open(out.c_str());
    
    if( !outputFile.is_open() ){
        std::cerr << "ERROR (write_to_file): unable to open for output file " << out << std::endl;
        exit(1);
    }
    
    if( rooted ){
        // must have odd number of sequences
        if( matrix.size() % 2 == 0  ){
            std::cerr << "ERROR (GHcaller): even number of sequences in rooted fasta matrix is wrong" << std::endl;
            exit(1);
        }
        
        for( unsigned int i = 0; i < matrix.size() - 1 ; i+=2){
            outputFile << ">" << names[i] << "\n";
            for ( int site = start0; site <= end0; site++) {
                outputFile << toIUPAC2( matrix[i][site] , matrix[i+1][site] ); //  matrix[i][site];
            }
            outputFile << std::endl;
        }
        outputFile << ">" << names[matrix.size() - 1] << "\n";
        for ( int site = start0; site <= end0; site++) {
            outputFile << toupper(matrix[matrix.size() - 1][site],loc); //  matrix[i][site];
        }
        outputFile << std::endl;
        
        
    }
    else{
        // must have even number of sequences
        if( matrix.size() % 2 != 0  ){
            std::cerr << "ERROR (GHcaller): odd number of sequences in unrooted fasta matrix is wrong" << std::endl;
            exit(1);
        }
        for( unsigned int i = 0; i < matrix.size() - 1 ; i+=2){
            outputFile << ">" << names[i] << "\n";
            for ( int site = start0; site <= end0; site++) {
                outputFile << toIUPAC2( matrix[i][site] , matrix[i+1][site] ); //  matrix[i][site];
            }
            outputFile << std::endl;
        }
        
        
    }
    
    outputFile.close();
    
}

// ADD REFERENCE/OUTGROUP

void fasta::append_seq( fasta in  ){
    // assuming fasta in to be reference, must be longer than fasta object
    // also must be a single sequence
    if( in.num_bases() < matrix.size() ){
        std::cerr << "ERROR (GHcaller): reference/outgroup used is shorter than data from mpileup ("
        << in.num_bases() << " and " << this->num_bases() << ")" << std::endl;
        exit(1);
    }
    else if ( in.num_bases() == this->num_bases()  ){
        matrix.push_back(in.matrix.at(0).c_str());
        names.push_back( in.name_at(0) );
    }
    else{
        std::string missing ( in.num_bases() - this->num_bases(), 'n');
        for (unsigned int nline = 0 ; nline < matrix.size() ; nline++) {
           matrix.at(nline).append(missing);
        }
        matrix.push_back(in.matrix.at(0).c_str());
        names.push_back( in.name_at(0) );
    }
    
}

// BUILD FROM GHcaller

void fasta::build_from_ngh( unsigned int site, std::string col_to_add ){
    
    while (matrix.at(0).size() < ( site - 1 ) ) {
        for (unsigned int i = 0; i < col_to_add.size(); i++) {
            matrix.at(i).push_back('n');
        }
    }
    
    for (unsigned int i = 0; i < col_to_add.size(); i++) {
        matrix.at(i).push_back(col_to_add.at(i));
    }
}

void fasta::set_all_names(std::vector<std::string> inames){
    // function should be called when matrix still has 2 lines per individual, before adding outgroup
    if( matrix.size() % 2 != 0  ){
        std::cerr << "ERROR (GHcaller): odd number of sequences in unrooted fasta matrix is wrong" << std::endl;
        exit(1);
    }
    
    std::stringstream aname;
    if( inames.size() == 1 &&  inames.at(0) == "NA" ){
        //std::cout << "adding default names to fasta file" << std::endl;
        for (unsigned int i = 0; i < matrix.size() / 2; i++) {
            aname << "i" << i + 1;
            names.push_back(aname.str());
            names.push_back(aname.str());
            aname.str("");
            
        }
        
    }
    else if ( inames.size() != matrix.size()/2 ){
        std::cerr << "WARNING (GHcaller): number of names provided (" << inames.size()
        << ") and individuals in input (" << matrix.size()/2
        << ") does not match. Will use default names" << std::endl;
        for (unsigned int i = 0; i < matrix.size() / 2; i++) {
            aname << "i" << i + 1;
            names.push_back(aname.str());
            names.push_back(aname.str());
            aname.str("");
            
        }
        
        
    }
    else {
        // correct number of names
        for (unsigned int n = 0; n < inames.size(); n++) {
            names.push_back(inames.at(n));
            names.push_back(inames.at(n));
        }
    }
    
    
}

void fasta::check_bed (std::string file_name){
    fasta_bed bed_regions;
    bed_regions.infile = file_name;
    std::string line;
    std::ifstream infile;
    std::stringstream ss;
    
    infile.open(file_name.c_str());
    if( !infile.is_open() ) {
        std::cerr << "ERROR (check_bed): Unable to open for reading bed file " << file_name << std::endl;
        exit(1);
    }
    
    while ( ! infile.eof() )
    {
        getline(infile, line);
        
        if (line.length() == 0)
            continue;
        std::vector <std::string> fields;
        msplit2( line, "\t", &fields );
        switch (fields.size()) {
            case 3:
                if ( bed_regions.chr.length() > 0 && bed_regions.chr != fields.at(0)){
                    std::cerr << "ERROR (check_bed): Bed file contains multiple chromosome, only 1 chromosome allowed! (" << bed_regions.chr << " and " << fields.at(0)  << ")"<<  std::endl;
                    exit(1);
  
                }
                else{
                    bed_regions.chr  = fields.at(0);
                }
                bed_regions.rstarts.push_back( atoi(fields.at(1).c_str()) );
                bed_regions.rends.push_back( atoi(fields.at(2).c_str()) );
                ss << "region_" << bed_regions.rstarts.size();
                bed_regions.rnames.push_back(ss.str() );
                ss.str("");
                break;
            case 2:
                std::cerr << "ERROR (check_bed): Malformatted region (missing fields) in bed file " << file_name << ", offending line: " << line <<  std::endl;
                exit(1);
                break;
            default:
                if ( bed_regions.chr.length() > 0 && bed_regions.chr != fields.at(0)){
                    std::cerr << "ERROR (check_bed): Bed file contains multiple chromosome, only 1 chromosome allowed! (" << bed_regions.chr << " and " << fields.at(0)  << ")"<<  std::endl;
                    exit(1);
                    
                }
                else{
                    bed_regions.chr  = fields.at(0);
                }
                // ignore optional columns in bed file
                bed_regions.rstarts.push_back( atoi(fields.at(1).c_str()) );
                bed_regions.rends.push_back( atoi(fields.at(2).c_str()) );
                bed_regions.rnames.push_back(fields.at(3));
                break;
        }
        
        if( bed_regions.rstarts.at(bed_regions.rstarts.size() -1 ) == 0
           || bed_regions.rends.at(bed_regions.rstarts.size() -1 ) == 0
           || bed_regions.rstarts.at(bed_regions.rstarts.size() -1 ) >= bed_regions.rends.at(bed_regions.rstarts.size() -1 )
           ){
            std::cerr << "ERROR (check_bed): Malformatted region in bed file " << file_name << ", offending line: " << line <<  std::endl;
            exit(1);
        }
        
        
    }
    infile.close();
    
    regions = bed_regions;
    
}



