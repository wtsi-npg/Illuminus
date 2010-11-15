////////////////////////////////////////////////////////////////////////
// ILLUMINUS version 2: a program that calls genotypes from the Illumina platform
// Copyright (c) 2007 Genome Research Ltd and University of Oxford
// Authors: Taane Clark <tc5@sanger.ac.uk, tgc@well.ox.ac.uk>, YY Teo <yy.teo@well.ox.ac.uk> 
// We would like to acknowledge Magnus Manske for his assistance in optimising aspects of the code, 
// as well as Mike Inouye for testing earlier versions. An earlier version of Illuminus used 
// some C++ code co-developed by Jonathan Marchini as a reference for data and variable
// input. This code was not in any way part of the key algorithmic development at the core of our 
// Illuminus software. The current version of Illuminus uses two non-standard C++ libraries
// (a) newmat11 (Robert Davies: http://www.robertnz.net/) and 
// (b) lapack++ (http://lapackpp.sourceforge.net/), both of which are free for commercial use. 
// Please report any problems to the authors above.
// Description: The code reads in a text file (columns: rs, coord, alleles, id_1a,
// id_1b, id_2a, id_2b, etc), and iterates using an EM algorithm to a 
// convergent set of calls. 
// There are several options: 
// -a perturbation analysis  
// -w whole genome amplified data
////////////////////////////////////////////////////////////////////////

//
// $Id$
//

using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string> 
#include <iomanip> 
#include <vector> 
#include <set> 
#include <map> 
#include <iterator> 
#include <algorithm> 
#include <bitset> 
#include <numeric> 
#include <functional> 
#include <cmath> 
#include <math.h>
#include "time.h"
#include "illuminus.h"

int version = 201;

int main(int argc, char *argv[]) {

  int i, seed1, seed2, rs = 1, niter = 100, low=-1 , upp= -1;
  char *infile = NULL, *outfile = NULL, *xfile = NULL, *yfile = NULL, *mfile = NULL;
  double thres = 0.95000001; 
  bool wga = false, pert = false, bed = false, calls = false, probs = false, pr = false;
  data dat;
  
  srand(time(NULL));
  seed1 = rand() % (10000) + 1;
  seed2 = rand() % (10000) + 1;
  
  //////////////////////
  // Input parameters //
  //////////////////////
  
  for(i = 0; i < argc; i++){
    if(*argv[i] == '-'){ 
			if ( strcmp ( argv[i] , "-nrs" ) == 0 ) { // don't use random seed
				seed1 = 12345;
				seed2 = 54321;
				continue ;
			}
      
      switch ( *(argv[i]+1) ) {
				case 'i' : infile = argv[++i]; break ; // data file 
				case 'x' : xfile = argv[++i]; break ; // chromosome X file
				case 'y' : yfile = argv[++i]; break ; // chromosome Y file .... not implemented
				case 'm' : mfile = argv[++i]; break ; // Mitochondrial file ..... not implemented
				case 'o' : outfile = argv[++i]; break ; // output file
				case 'n' : niter = atoi(argv[++i]); break ; // max number of iterations 
				case 't' : thres = atof(argv[++i]); break ; // trimming threshold for calling probabilities 
				case 'c' : calls = true; break ; // write calls to file outfile_c
				case 'p' : probs = true; break ; // write probs of calls to outfile_p
				case 'b' : bed = true; break ; // write calls to plink bed
				case 'w' : wga = true; break ; // wga
				case 'a' : pert = true;  break ; // perturbation analysis
				case 'v' : cout << endl << "Illuminus version " << version/100.0 << endl 
				                << "Compiled: " << __DATE__ << endl; 
				           exit(0);
				case 's' : { 
					low = atoi(argv[++i]);
					upp = atoi(argv[++i]); /// range of snps
					break ;
				}
      }
    }
  }

  if( outfile == NULL or (infile == NULL and xfile == NULL and yfile == NULL and mfile == NULL)) { 
    cout << "Need to specify both an input and an output file" << endl; exit(1); }        
  if( calls == false and probs == false) { cout << "Need to specify -c and / or -p" << endl; exit(1); }     
  
  /////////////////////////////////
  // Set random number generator //
  /////////////////////////////////
  
  setall(seed1, seed2);
  
  ////////////////////////////////////////
  // Inputs to the class and ruuning it //
  ////////////////////////////////////////
  
  dat = data(infile, xfile, low, upp);
  dat.outfile = outfile;
  dat.niter = niter;
  dat.thres = thres;
  dat.niter = niter;
  dat.calls = calls;
  dat.bed = bed;
  dat.probs = probs;
  dat.wga = wga;
  dat.pert = pert;
 
  dat.illuminus();
  
}

//// example of use:   ./illuminus -i example.txt -o out -c -a -p 
//// this will run illuminus on example.txt, outputting out_calls and out_probs,
//// as well as performing a perturbation analysis reported in both files
