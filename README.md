Genotype/Haplotype SNP caller 
 
Authors: B. Nevado & S. Ramos-Onsins (used in B. Nevado, S.E. Ramos-Onsins and M. Perez-Enciso Mol.Ecol. 2014) 
 
Installation:  
 
Linux: g++ -Wall -g -std=c++0x -O3 *.cpp -o GHcaller  
Mac: clang++ -std=c++11 -stdlib=libc++ -O3 -Wall *.cpp -o GHcaller (for Mac, clang++ needs to be version 4+) 
 
 
Description: 
 
Genotype-Haplotype SNP calling from mpileup data  
Input: mpileup data from STDIN, for a single chromosome, 
  e.g. "samtools mpileup -r chr1 infile1.bam infile2.bam > data.chr1.mpileup" 
Output: fasta file 
 
 
Usage: cat data.mpileup | GHcaller -outfile outfile [ filtering options ] [ SNP calling options ] [ Output options ] 
 
[ filtering options ]  
  -baseq: minimum base quality of reads (def: 20). Reads below threshold are discarded.  
  -mindep: minimum number of reads (def: 3). Sites below threshold are coded missing.  
  -maxdep: maximum number of reads (def: 100). Sites above threshold are coded missing.  
  -platform: offset to convert base qualities, e.g. 33 is for illumina 1.8+ and Sanger, 64 for Illumina 1.3, ... (def: 33)  
  *Filtering options accept comma-separated list of values per individual, or single value for all individuals.  
   
[ SNP calling options ]  
  -minreads: minimum number of reads to attempt a genotype call (def: 6).  
  -haplotypes: set to zero (0) to perform only genotype calls (def: 1, i.e. perform both haplotype and genotype calls).  
  -error: raw sequencing error (def: 0.01).  
  -chivalue: chi-square threshold for LRT test of genotypes (def: 3.841459, i.e. 0.05 with 1 d.f.).  
 
[ Output options ]  
  -verbose: set to 0 to supress info output to screen (def: 1).  
  -outfile: file to write to (required).  
  -names: comma-separated list of individuals' names, in same order as input (def: individuals are named i1..i2..in).  
  -outgroup|-reference: sequence file in fasta format to add to output alignment (def: none).  
  -iupac: if 1, will output 1 sequence per individual using IUPAC codes, and lower-case bases denoting haplotype calls (def: 0, output 2 lines per individual)  
  -bed: output regions in bed file only (def: none)  
  -concatenate: set to 1 to concatenate all regions in input bed into a single output file (def: 0, output 1 file per region in bed file)  
  -wlen: output fasta sequences in windows of specified length, 1 window per file (def: none) 
 
Algorithm is based on :  
Lynch (2009) Genetics 182: 295-301  
Roesti et al. (2012) Molecular Ecology 21: 2852–2862. 
 
This code has been partially developed thanks to the Grant CGL2009-09346 (MICINN, Spain).  
Also available from : http://bioinformatics.cragenomica.es/numgenomics/people/sebas/software/software.html  
Contact: bruno.nevado@gmail.com 
