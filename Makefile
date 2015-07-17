
PROGRAM=GHcaller

CC=g++
CFLAGS=-Wall -O3 -std=c++11

$PROGRAM: source/common.cpp  source/common.h  source/fasta.cpp  source/fasta.h  source/main.cpp  source/ngh.cpp  source/ngh.h
	$(CC) $(CFLAGS) -o $(PROGRAM)  source/common.cpp  source/common.h  source/fasta.cpp  source/fasta.h  source/main.cpp  source/ngh.cpp  source/ngh.h
