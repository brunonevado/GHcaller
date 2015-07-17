//
//  ngh.cpp
//  GHcaller
//
//  Copyright (c) 2013 Bruno Nevado, GNU license.

#include <cctype>
#include <set>
#include <string>
#include <cstring>
#include <cmath>  // lgamma, log
#include <algorithm>
#include <locale>
#include <sstream>

#include "ngh.h"

ngh::ngh( std::string mindep, std::string maxdep, std::string platform, std::string baseq,
         int mreads, float err, float chi, int wlen, int iiupac, bool haplos ){
    DNAbases = { 'a', 'A', 'c', 'C', 'g', 'G', 't', 'T' };
    // AA, AC, AG, AT, CC, CG, CT, GG, GT, TT, for each ind
    
    std::vector <std::string> temp;
    msplit2(mindep, ",", & temp);
    for(unsigned int i = 0; i < temp.size(); i++)
        mindeps.push_back( atoi ( temp.at(i).c_str() ) );
    temp.clear();
    
    msplit2(maxdep, ",", &temp);
    for(unsigned int i = 0; i < temp.size(); i++)
        maxdeps.push_back( atoi ( temp.at(i).c_str() ) );
    temp.clear();
    
    msplit2(platform, ",", & temp);
    for(unsigned int i = 0; i < temp.size(); i++)
        platforms.push_back( atoi ( temp.at(i).c_str() ) );
    temp.clear();
    
    msplit2(baseq, ",", & temp);
    for(unsigned int i = 0; i < temp.size(); i++)
        baseqs.push_back( atoi ( temp.at(i).c_str() ) );
    temp.clear();
    
    minreads = mreads;
    p_err = err;
    chivalue = chi;
    winlen = wlen;
    iupac = iiupac;
    do_haplos = haplos;
}

void ngh::set_num_inds (  unsigned int n ){
    num_inds = n;
    
    
    
    //  min deps
    if( ( mindeps.size() != 1 && n > 1 )  && mindeps.size() != n ){
        std::cerr << "ERROR (ngh): number of mindeps provided (" << mindeps.size()
        << ") and individuals in input (" << n << ") does not match\n";
        exit(1);
    }
    else if ( mindeps.size() == 1  && n > 1){
        mindeps.resize(n, mindeps.at(0) );
        
    }
    
    //  max deps
    if( ( maxdeps.size() != 1 && n > 1 )  && maxdeps.size() != n ){
        std::cerr << "ERROR (ngh): number of maxdeps provided (" << maxdeps.size()
        << ") and individuals in input (" << n << ") does not match\n";
        exit(1);
    }
    else if ( maxdeps.size() == 1  && n > 1){
        maxdeps.resize(n,maxdeps.at(0));
    }
    
    //  platforms
    if( ( platforms.size() != 1 && n > 1 )  && platforms.size() != n ){
        std::cerr << "ERROR (ngh): number of platforms provided (" << platforms.size()
        << ") and individuals in input (" << n << ") does not match\n";
        exit(1);
    }
    else if ( platforms.size() == 1  && n > 1){
        platforms.resize(n,platforms.at(0));
        
        
    }
    
    //  min base qualities
    if( ( baseqs.size() != 1 && n > 1 )  && baseqs.size() != n ){
        std::cerr << "ERROR (ngh): number of baseqs provided (" << baseqs.size()
        << ") and individuals in input (" << n << ") does not match\n";
        exit(1);
    }
    else if ( baseqs.size() == 1  && n > 1){
        baseqs.resize(n,baseqs.at(0));
        
    }
    
}

std::vector <std::string> ngh::info_to_vector(){
    
    std::stringstream filtering, snpcallinfo;
    filtering << "Num inds: " << num_inds <<  ", mindeps:";
    for(unsigned int i = 0; i < mindeps.size(); i++)
        filtering << " " << mindeps.at(i);
    
    filtering << ", maxdeps:";
    for(unsigned int i = 0; i < maxdeps.size(); i++)
        filtering << " " << maxdeps.at(i);
    
    filtering << ", platforms:";
    for(unsigned int i = 0; i < platforms.size(); i++)
        filtering  << " " << platforms.at(i);
    
    filtering << ", baseqs:";
    for(unsigned int i = 0; i < baseqs.size(); i++)
        filtering << " " << baseqs.at(i);
    
    snpcallinfo << "SNP calling algorithm: Roesti, minreads for genotype: " << minreads
      << ", error: " << p_err << ", LRT chi-square threshold: " << chivalue;
    if(do_haplos)
        snpcallinfo << ", haplotype calls: enabled";
    else
        snpcallinfo << ", haplotype calls: disabled";
            
    
    std::vector <std::string > toreturn;
    toreturn.push_back(filtering.str());
    toreturn.push_back(snpcallinfo.str());
    return toreturn;
}

ngh_parse_return ngh::parse_mpileup( std::string line ){
    
    std::locale loc;
    
    unsigned int site;
    std::string chr;
    char ref;
    std::vector <std::string> vreads, vquals;
    std::vector <int> raw_deps;
    std::vector <std::string> fields;
    msplit2(line, "\t" , &fields);
    
    unsigned int ninds = ( (fields.size() - 3 ) / 3 );
    site = atoi ( fields.at(1).c_str() ) ;
    chr = fields.at(0);
    ref = tolower( fields.at(2).at(0), loc) ;
    
    for ( unsigned int ind0 = 0; ind0 < ninds ; ind0++) {
        unsigned int ind0depth = atoi ( fields.at( ind0 * 2 + ind0 + 3 ).c_str() );
        raw_deps.push_back(  ind0depth );
        std::string reads = fields.at( ind0 * 2 + ind0 + 4 );
        std::string quals = fields.at( ind0 * 2 + ind0 + 5 );
        
        // 1. remove beginning/end of reads and indels info
        for( int del = ( int (reads.size()) - 1 ) ; del >= 0 ; del-- ){
            if ( reads.at(del) == '^' ){
                // means beginning of read, next symbol is mapping quality in ascii - remove both
                if( del > 0 && reads.at(del - 1 ) == '^' )
                    continue;
                
                reads.erase( reads.begin()+del );
                reads.erase( reads.begin()+del, reads.begin()+del+1 );
                continue;
            }
        }
        for( int del = ( int (reads.size()) - 1 ) ; del >= 0 ; del-- ){
            
            if( reads.at(del) == '.' ){
                reads.at(del) = toupper(ref);
                continue;
            }
            else if ( reads.at(del) == ','){
                reads.at(del) = tolower(ref);
                continue;
            }
            else if ( reads.at(del) == '$' ){
                reads.erase( reads.begin()+del );
                continue;
            }
            
            else if ( reads.at(del) == '+' || reads.at(del) == '-'  ){
                unsigned int indel_size = abs( atoi ( reads.substr( del ).c_str()  ));
                unsigned int num_length =  int (strlen( std::to_string(  static_cast<long long> (indel_size) ).c_str() ));
                reads.erase( reads.begin()+del, reads.begin()+del+1+num_length + indel_size );
                continue;
                
            }
            
        }
        
        // 2. filter by quality (also remove bases marked '*' which mean insertion)
        for( int del = ( int (reads.size()) - 1 ) ; del >= 0 ; del-- ){
            if( reads.at(del) == '*' || (((int) quals.at(del)) - platforms.at(ind0) ) < baseqs.at(ind0)  ){
                reads.erase(reads.begin()+del);
                quals.erase(quals.begin()+del);
                
            }
            
        }
        if ( reads.size() < mindeps.at(ind0) || reads.size() > maxdeps.at(ind0)){
            vreads.push_back("N");
            vquals.push_back("N");
        }
        else{
            vreads.push_back(reads);
            vquals.push_back(quals);
            
            // this check is a pain?
            if( reads.size() != quals.size() ){
                std::cerr << "ERROR (parse mpileup): reads and quals dont match after parsing, site " << site
                << " ind " << ind0+1 << " reads: " << reads << " quals: " << quals << std::endl;
                exit(1);
            
            }
        }
    
    
    
    
    
    
    }
    
    
    ngh_parse_return toreturn;
    toreturn.chr = chr;
    toreturn.site = site;
    toreturn.ninds = ninds;
    toreturn.ref = ref;
    toreturn.quals = vquals;
    toreturn.reads = vreads;
    toreturn.raw_deps = raw_deps;
     
    return toreturn;
    
}

res_snp_call ngh::roesti_snpcall ( ngh_parse_return in ){
    
    res_snp_call toreturn;
    std::vector < std::vector <float>> gls;
    for (unsigned int ind0 = 0; ind0 < num_inds; ind0++ ) {
        if( in.reads.at(ind0).size() < 2 ){
            toreturn.called_genotypes.push_back('n');
            toreturn.called_genotypes.push_back('n');
            std::vector <float> ps(10, 0.0);
            gls.push_back(ps);
            continue;
        }
        else if ( in.reads.at(ind0).size() < minreads  ){
            // will do haplocall
            std::vector <float> ps = calc_genos_probs ( in.reads.at(ind0));
            gls.push_back(ps);
            if( !do_haplos ){
                toreturn.called_genotypes.push_back('n');
                toreturn.called_genotypes.push_back('n');
                continue;
                
            }
            int best_haplo = lrt_haplotypes(ps);
            switch ( best_haplo ) {
                case -1:
                    toreturn.called_genotypes.push_back( 'n' );
                    toreturn.called_genotypes.push_back( 'n' );
                    break;
                case 0:
                    toreturn.called_genotypes.push_back( 'a' );
                    toreturn.called_genotypes.push_back( 'n' );
                    break;
                case 1:
                    toreturn.called_genotypes.push_back( 'c' );
                    toreturn.called_genotypes.push_back( 'n' );
                    break;
                case 2:
                    toreturn.called_genotypes.push_back( 'g' );
                    toreturn.called_genotypes.push_back( 'n' );
                    break;
                case 3:
                    toreturn.called_genotypes.push_back( 't' );
                    toreturn.called_genotypes.push_back( 'n' );
                    break;
                    
                default:
                    std::cerr << "ERROR (GHcaller): error in return from lrt_haplotypes" << std::endl;
                    exit(1);
                    break;
            }
            
        }
        else{
            // attempt genotype call
            std::vector <float> ps = calc_genos_probs ( in.reads.at(ind0));
            gls.push_back(ps);
            
            int best_genotype = lrt_genotypes ( ps );
            if( !do_haplos && best_genotype == -1 ){
                toreturn.called_genotypes.push_back('n');
                toreturn.called_genotypes.push_back('n');
                continue;
            }
            switch ( best_genotype ) {
                case -1:
                    switch ( lrt_haplotypes(ps) ) {
                        case -1:
                            toreturn.called_genotypes.push_back( 'n' );
                            toreturn.called_genotypes.push_back( 'n' );
                            break;
                        case 0:
                            toreturn.called_genotypes.push_back( 'a' );
                            toreturn.called_genotypes.push_back( 'n' );
                            break;
                        case 1:
                            toreturn.called_genotypes.push_back( 'c' );
                            toreturn.called_genotypes.push_back( 'n' );
                            break;
                        case 2:
                            toreturn.called_genotypes.push_back( 'g' );
                            toreturn.called_genotypes.push_back( 'n' );
                            break;
                        case 3:
                            toreturn.called_genotypes.push_back( 't' );
                            toreturn.called_genotypes.push_back( 'n' );
                            break;
                        default:
                            std::cerr << "ERROR (GHcaller): error in return from lrt_haplotypes" << std::endl;
                            exit(1);
                            break;
                    }
                    break;
                case 0:
                    toreturn.called_genotypes.push_back( 'a' );
                    toreturn.called_genotypes.push_back( 'a' );
                    break;
                case 1:
                    toreturn.called_genotypes.push_back( 'a' );
                    toreturn.called_genotypes.push_back( 'c' );
                    break;
                case 2:
                    toreturn.called_genotypes.push_back( 'a' );
                    toreturn.called_genotypes.push_back( 'g' );
                    break;
                case 3:
                    toreturn.called_genotypes.push_back( 'a' );
                    toreturn.called_genotypes.push_back( 't' );
                    break;
                case 4:
                    toreturn.called_genotypes.push_back( 'c' );
                    toreturn.called_genotypes.push_back( 'c' );
                    break;
                case 5:
                    toreturn.called_genotypes.push_back( 'c' );
                    toreturn.called_genotypes.push_back( 'g' );
                    break;
                case 6:
                    toreturn.called_genotypes.push_back( 'c' );
                    toreturn.called_genotypes.push_back( 't' );
                    break;
                case 7:
                    toreturn.called_genotypes.push_back( 'g' );
                    toreturn.called_genotypes.push_back( 'g' );
                    break;
                case 8:
                    toreturn.called_genotypes.push_back( 'g' );
                    toreturn.called_genotypes.push_back( 't' );
                    break;
                case 9:
                    toreturn.called_genotypes.push_back( 't' );
                    toreturn.called_genotypes.push_back( 't' );
                    break;
                default:
                    std::cerr << "ERROR (GHcaller): error in return from lrt_genotypes" << std::endl;
                    exit(1);
                    
                    break;
            }
            
        }
        
        
    }
    toreturn.gls = gls;
    return toreturn;
    
}

int ngh::lrt_genotypes ( std::vector <float> in  ) {
    // returns index of most likely genotype (in AA, AC, AG, AT, CC, CG, CT, GG, GT, TT)
    // or -1 if LRT between 2 most likely genotypes is not significant
    float maxGL = *max_element (in.begin(),in.end());
    
    int index_max = -1;
    for (unsigned int i = 0; i < in.size(); i++) {
        if (in.at(i) == maxGL){
            in.erase( in.begin() + i);
            index_max = i;
            break;
        }
    }
    
    float smaxGL = *max_element (in.begin(),in.end());
    float lrt = abs( 2 * ( maxGL - smaxGL) );
    if(lrt < chivalue){
        return -1;
        
    }
    else{
        return index_max;
    }
}

int ngh::lrt_haplotypes ( std::vector <float> in_all  ) {
    // returns index of most likely haplotype (in AA, CC, GG, TT)
    // or -1 if LRT between 2 most likely haplotypes is not significant
    std::vector < float > in { in_all.at(0), in_all.at(4), in_all.at(7), in_all.at(9) };
    float maxGL = *max_element (in.begin(),in.end());
    int index_max = -1;
    
    for (unsigned int i = 0; i < in.size(); i++) {
        if (in.at(i) == maxGL){
            in.erase( in.begin() + i);
            index_max = i;
            break;
        }
    }
    
    float smaxGL = *max_element (in.begin(),in.end());
    float lrt = abs( 2 * ( maxGL - smaxGL) );
    if(lrt < chivalue){
        return -1;
    }
    else{
    return index_max;
    }
}


std::vector <float> ngh::calc_genos_probs ( std::string reads ){
    
    std::vector < float > toreturn;
    // calculates GLs for AA, AC, AG, AT, CC, CG, CT, GG, GT, TT
    // DNAbases = { 'a', 'A', 'c', 'C', 'g', 'G', 't', 'T' };
    
    // count each base in forward and reverse strands
    std::vector <int> counts = { 0, 0, 0, 0, 0, 0, 0, 0 };
    for( unsigned int each_read = 0; each_read < reads.size();each_read++ ){
        for( unsigned int each_base = 0; each_base < DNAbases.size(); each_base++ ){
            if( reads.at(each_read) == DNAbases.at(each_base) ){
                counts.at(each_base)++;
            }
            
        }
    }
    
    int fA = counts.at(0) + counts.at(1);
    int fC = counts.at(2) + counts.at(3);
    int fG = counts.at(4) + counts.at(5);
    int fT = counts.at(6) + counts.at(7);
    
    toreturn.push_back(calc_homo_p( fA, fC, fG, fT));
    toreturn.push_back(calc_hetero_p( fA, fC, fG, fT));
    toreturn.push_back(calc_hetero_p( fA, fG, fC, fT));
    toreturn.push_back(calc_hetero_p( fA, fT, fC, fG));
    toreturn.push_back(calc_homo_p( fC, fA, fG, fT));
    toreturn.push_back(calc_hetero_p( fC, fG, fA, fT));
    toreturn.push_back(calc_hetero_p( fC, fT, fG, fA));
    toreturn.push_back(calc_homo_p( fG, fA, fC, fT));
    toreturn.push_back(calc_hetero_p( fG, fT, fC, fA));
    toreturn.push_back(calc_homo_p( fT, fA, fG, fC));
    
    return toreturn;
}


float ngh::calc_homo_p ( int n1, int n2, int n3, int n4 ){
    
    float lg1 = lgamma( n1 + 1 ) ;
    float lg2 = lgamma( n2 + 1 ) ;
    float lg3 = lgamma( n3 + 1 ) ;
    float lg4 = lgamma( n4 + 1 ) ;
    float lgt = lgamma( n1 + n2 + n3 + n4 + 1 );
    float res = (lgt-lg1-lg2-lg3-lg4+(n1)*log(1.0-3.0/4.0*p_err)+( n2 + n3 + n4 )*log( p_err /4.0));
    return res;
    
}


float ngh::calc_hetero_p ( int n1, int n2, int n3, int n4 ){
    
    float lg1 = lgamma( n1 + 1 );
    float lg2 = lgamma( n2 + 1 );
    float lg3 = lgamma( n3 + 1 );
    float lg4 = lgamma( n4 + 1 );
    float lgt = lgamma( n1 + n2 + n3 + n4 + 1 );
    float res = (lgt-lg1-lg2-lg3-lg4+(n1+n2)*log(0.5/2.0-p_err/4.0)+(n3+n4)*log(p_err/4.0));
    return res;
    
}
