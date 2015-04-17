//
//  main.cpp
//  GHcaller
//
//  Copyright (c) 2013 Bruno Nevado, GNU license.

//  for linux, compile with g++ -Wall -g -std=c++0x -O3  *.cpp -o GHcaller
//  for mac with clang++ -std=c++11 -stdlib=libc++ -O3 -Wall *.cpp -o GHcaller

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>

#include "ngh.h"
#include "fasta.h"


int main(int argc, char* argv[])
{
    
    
    std::string version = "<GHcaller> v 0.0.1 05122013";
    std::string line;
    std::vector <std::string> args;
    std::vector <std::string> values;
    
    if( argc == 2 && !strncmp (argv[1], "-help", 5) ){
        help(version);
        exit(0);
    }
    if ( argc % 2 == 0 ){
        
        std::cerr << "ERROR(GHcaller): wrong argument number" << std::endl;
        help( version );
        exit(1);
    }
    for( int i = 1; i < argc; i += 2 )
    {
        args.push_back(argv[i]);
        values.push_back(argv[i+1]);
    }
    
    // FILTERING OPTIONS
    std::string mindep = "3",
    maxdep = "100",
    platform = "33",
    baseq = "20";
    
    // SNP CALLING ALGORITHM OPTIONS
    int minreads = 6;  // genotype / haplotype call threshold
    float p_err = 0.01,   // raw sequencing error
    chivalue = 3.841459 ; // 1 d.f., 0.05 (for LRT)
    bool haplocalls = true;  // to disable haplotpye calls
    
    // OUTPUT
    int iupac = 0,      // if 1 uses iupac + upper/lower case bases
    wlen = 0;           // window len
    std::string outgroup_file = "NA",  // will append sequence to outfile (and fill ends of seqs with Ns)
    outfile = "NA",
    bed_file = "NA";
    bool verbose = true,
    concatenate_fastas = false;   // when using bed file only
    
    std::string allnames = "NA";
    
    // GET ARGS
    for(unsigned int i = 0; i < args.size(); i++ ){
        if( strncmp (args.at(i).c_str(), "-", 1) ){
            std::cerr << "ERROR (GHcaller): fail to understand argument " << args.at(i) << "\n";
            exit(1);
        }
        // FILTERING OPTIONS
        else if( !strncmp (args.at(i).c_str(), "-mindep", 5)  ){
            mindep = values.at(i) ;
        }
        else if( !strncmp (args.at(i).c_str(), "-maxdep", 5)  ){
            maxdep = values.at(i) ;
        }
        else if( !strncmp (args.at(i).c_str(), "-platform", 4)  ){
            platform = values.at(i) ;
        }
        else if( !strncmp (args.at(i).c_str(), "-baseq", 4)  ){
            baseq = values.at(i) ;
        }
        
        // SNP CALL OPTIONS
        else if( !strncmp (args.at(i).c_str(), "-haplotypes", 5) ){
            if( values.at(i) == "1"  ){
                haplocalls = true;
            }
            else if(  values.at(i) == "0"  ){
                haplocalls = false;
            }
            else{
                std::cerr << "ERROR (GHcaller): unknown haplotype calls option (" << values.at(i)
                << "). Available values are 1 (do haplotype calls) and 0 (skip haplotype calls)" << std::endl;
                exit(1);
            }
        }
        else if( !strncmp (args.at(i).c_str(), "-minreads", 5) ){
            minreads = atoi( values.at(i).c_str() ) ;
        }
        else if(!strncmp (args.at(i).c_str(), "-error", 3) ){
            p_err = std::stof(values.at(i).c_str());
        }
        else if(!strncmp (args.at(i).c_str(), "-chivalue", 3) ){
            chivalue = std::stof(values.at(i).c_str());
        }
        // OUTPUT OPTIONS
        else if( !strncmp (args.at(i).c_str(), "-outfile", 5) ){
            outfile = values.at(i) ;
        }
        else if( !strncmp (args.at(i).c_str(), "-outgroup", 5) || !strncmp (args.at(i).c_str(), "-reference", 4) ){
            outgroup_file = values.at(i) ;
        }
        else if( !strncmp (args.at(i).c_str(), "-bed", 5) ){
            bed_file = values.at(i) ;
        }
        else if( !strncmp (args.at(i).c_str(), "-winlen", 5) ){
            wlen = atoi( values.at(i).c_str() ) ;
        }
        else if(!strncmp (args.at(i).c_str(), "-concatenate", 3) ){
            concatenate_fastas = std::stoi(values.at(i).c_str());
        }
        else if(!strncmp (args.at(i).c_str(), "-iupac", 3) ){
            iupac = std::stoi(values.at(i).c_str());
        }
        else if(!strncmp (args.at(i).c_str(), "-names", 3) ){
            allnames = values.at(i);
        }
        else if(!strncmp (args.at(i).c_str(), "-verbose", 3) ){
            verbose = std::stoi(values.at(i).c_str());
        }
        else{
            std::cerr << "ERROR (GHcaller): fail to understand argument value pair: " << args.at(i) << " " << values.at(i) << std::endl;
            exit(1);
        }
        
    }
    
    // CHECK ARGS
    if(outfile == "NA" ){
        std::cerr << "ERROR (GHcaller): output file unspecified (-outfile)" << std::endl;
        exit(1);
    }
    
    if( wlen != 0 && bed_file != "NA" ){
        std::cerr << "ERROR (fasta output): -winlen and -bed options are mutually exclusive" << std::endl;
        exit(1);
    }
    if( concatenate_fastas == true && bed_file == "NA"  ){
        std::cerr << "ERROR (fasta output): -concatenate option only valid with bed file" << std::endl;
        exit(1);
    }
    
    std::vector <std::string> vnames;
    if( allnames != "NA"){
        msplit2(allnames, ",", & vnames);
    }
    else{
        vnames.push_back("NA");
    }
    
    if(verbose)
        std::clog << version << ", waiting for input from STDIN..." << std::endl;
    
    
    
    ngh angh( mindep, maxdep, platform, baseq, minreads, p_err, chivalue, wlen, iupac, haplocalls);
    ngh_parse_return parsed_data;
    
    std::string lineInput;
      
    getline(std::cin,lineInput);
    
    std::vector <std::string> fields;
    msplit2(lineInput, "\t" , &fields);
    
    
    angh.set_num_inds( (fields.size() - 3 )/3 );
    
    
    std::vector < std::string > info = angh.info_to_vector();
    
    if(verbose){
        for (unsigned int i = 0; i < info.size(); i++) {
            std::clog << "<GHcaller> " << info.at(i) << std::endl;
        }
        std::clog << "<GHcaller> Output format: fasta";
        if(iupac != 0)
            std::clog << ", using iupac codes";
        else
            std::clog << ", using unphased haplotypes";
        if (wlen != 0)
            std::clog << ", window size: " << wlen;
        if( bed_file != "NA")
            std::clog << ", bed file: " << bed_file;
        
        std::clog << ", output file: " << outfile;
        
        if(bed_file != "NA" && concatenate_fastas)
            std::clog << " (all regions)";
        else if( bed_file != "NA" && !concatenate_fastas)
            std::clog << " (prefix for each region)";
        
        std::clog  << std::endl;
    }
    
    
    fasta nameless_fasta(1);
    nameless_fasta.set_num_inds(2*(fields.size() - 3 )/3 );
    nameless_fasta.set_all_names(vnames);
    if( bed_file != "NA" )
        nameless_fasta.check_bed(bed_file);
 
    
    
    
    
    parsed_data = angh.parse_mpileup(lineInput);
    nameless_fasta.build_from_ngh( parsed_data.site,  angh.roesti_snpcall(parsed_data).called_genotypes );
    int nlines = 1;
    while (getline(std::cin,lineInput)) {
        nlines++;
        if( nlines % 10000000 == 0)
            std::clog << "<GHcaller> " << nlines << " lines processed" << std::endl;
        parsed_data = angh.parse_mpileup(lineInput);        
        nameless_fasta.build_from_ngh( parsed_data.site,  angh.roesti_snpcall(parsed_data).called_genotypes );
    }
    bool rooted = false;
    if( outgroup_file != "NA" ){
        fasta outgroup(1);
        outgroup.read_fasta_file(outgroup_file);
        if( outgroup.num_lines() > 1 ){
            std::cerr << "ERROR (GHcaller): outgroup file should contain a single sequence ("
            << outgroup.num_lines() << " found)" << std::endl;
            
        }
        else{
            nameless_fasta.append_seq( outgroup );
            if(verbose)
                std::clog << "<GHcaller> Rooted output with sequence from " << outgroup_file << std::endl;
        
            rooted = true;
        }
    }
    
    if( bed_file == "NA" )
        nameless_fasta.write(outfile, wlen, rooted, iupac, verbose );
    else
        nameless_fasta.write_bed(outfile, wlen, rooted, iupac, concatenate_fastas, verbose);
    
    
    return (0);
}
