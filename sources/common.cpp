//
//  common.cpp
//  NGHcaller
//
//  Copyright (c) 2013 Bruno Nevado, GNU license.

#include <fstream>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "common.h"

char toIUPAC2 (char in1, char in2 ){
    std::locale loc;
    
    // iupac 2 includes lower case for haplotype, and upper case for genotype
    if( in1 == in2 ){
        return toupper(in1,loc);
    }
    else if ( in1 == 'n' ){
        return ( in2 );
    }
    else if ( in2 == 'n' ){
        return ( in1 );
    }
    else if ( (in1 == 'a' && in2 == 'c') || (in2 == 'a' && in1 == 'c') ){
        return 'M';
    }
    else if ( (in1 == 'a' && in2 == 'g') || (in2 == 'a' && in1 == 'g') ){
        return 'R';
    }
    else if ( (in1 == 'a' && in2 == 't') || (in2 == 'a' && in1 == 't') ){
        return 'W';
    }
    else if ( (in1 == 'c' && in2 == 'g') || (in2 == 'c' && in1 == 'g') ){
        return 'S';
    }
    else if ( (in1 == 't' && in2 == 'g') || (in2 == 't' && in1 == 'g') ){
        return 'K';
    }
    else if ( (in1 == 'c' && in2 == 't') || (in2 == 'c' && in1 == 't') ){
        return 'Y';
    }
    else{
        std::cerr << "ERROR (turnIUPAC): no code available for " << in1 << in2 << "\n";
        exit(1);
    }
    
}

void msplit2( const std::string& s , std::string delim,  std::vector<std::string> * output){
    
    
    unsigned long start = 0U;
    unsigned long end = s.find(delim);
    while (end != std::string::npos)
    {
        output->push_back(s.substr(start, end - start));
        start = end + delim.length();
        end = s.find(delim, start);
    }
    output->push_back(s.substr(start, end - start));
    
}

void help( std::string v ){
    
    std::cout << std::endl << "Naive Genotype/Haplotype SNP caller"<< std::endl;
    std::cout << "Authors: B. Nevado & S. Ramos-Onsins"<< std::endl;
    std::cout << "Contact: bruno.nevado@cragenomica.es"<< std::endl;
    
    
    std::cout << "  Genotype-Haplotype SNP calling from mpileup data" << std::endl;
    
    
    std::cout << "  Input: mpileup data from STDIN, for a single chromosome, e.g. \"samtools mpileup -r chr1 infile1.bam infile2.bam > data.chr1.mpileup\"" << std::endl;
    std::cout << "  Output: fasta file" << std::endl << std::endl;
    
    std::cout << "  Usage: cat data.mpileup  | ./NGHcaller -outfile outfile [ filtering options ] [ SNP calling options ] [ Output options ]" << std::endl << std::endl;
    std::cout << "  [ filtering options ]" << std::endl;
    std::cout << "    -baseq: minimum base quality of reads (def: 20). Reads below threshold are discarded." << std::endl;
    std::cout << "    -mindep: minimum number of reads (def: 3). Sites below threshold are coded missing."   << std::endl;
    std::cout << "    -maxdep: maximum number of reads (def: 100). Sites above threshold are coded missing." << std::endl;
    std::cout << "    -platform: offset to convert base qualities, e.g. 33 is for illumina 1.8+ and Sanger, 64 for Illumina 1.3, ... (def: 33)" << std::endl;
    std::cout << "    *Filtering options accept comma-separated list of values per individual, or single value for all individuals." << std::endl;
    std::cout << "  [ SNP calling options ]" << std::endl;
    std::cout << "    -minreads: minimum number of reads to attempt a genotype call (def: 6)." << std::endl;
    std::cout << "    -haplotypes: set to zero (0) to perform only genotype calls (def: 1, i.e. perform both haplotype and genotype calls)." << std::endl;
    std::cout << "    -error: raw sequencing error (def: 0.01)." << std::endl;
    std::cout << "    -chivalue: chi-square threshold for LRT test of genotypes (def: 3.841459, i.e. 0.05 with 1 d.f.)." << std::endl;

    
    std::cout << "  [ Output options ]" << std::endl;
    
    std::cout << "    -verbose: set to 0 to supress info output to screen (def: 1)." << std::endl ;
    std::cout << "    -outfile: file to write to (required)." << std::endl;
    std::cout << "    -names: comma-separated list of individuals' names, in same order as input (def: individuals are named i1..i2..in)." << std::endl;
    std::cout << "    -outgroup|-reference: sequence file in fasta format to add to output alignment (def: none)." << std::endl;
    std::cout << "    -iupac: if 1, will output 1 sequence per individual using IUPAC codes, and lower-case bases denoting haplotype calls (def: 0, output 2 lines per individual)" << std::endl;
    std::cout << "    -bed: output regions in bed file only (def: none)" << std::endl;
    std::cout << "    -concatenate: set to 1 to concatenate all regions in input bed into a single output file (def: 0, output 1 file per region in bed file)" << std::endl;
    std::cout << "    -wlen: output fasta sequences in windows of specified length, 1 window per file (def: none)" << std::endl;
    
    std::cout  << std::endl << v << std::endl;
    
}

