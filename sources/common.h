//
//  common.h
//  NGHcaller
//
//  Copyright (c) 2013 Bruno Nevado, GNU license.


#ifndef __ngh_common__
#define __ngh_common__

#include <iostream>
#include <vector>
#include <string>

void msplit2( const std::string& tosplit , std::string delim,  std::vector<std::string> * output );

char  toIUPAC2 ( char in1, char in2 ); // for haplo/genotype calls, returns e.g. c for CN, and C for CC

void help ( std::string v );

#endif