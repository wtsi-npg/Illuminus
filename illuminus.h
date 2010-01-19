// ILLUMINUS version 2: a program that calls genotypes from the Illumina platform
// Copyright (c) 2007 Genome Research Ltd and University of Oxford
// Authors: Taane Clark <tc5@sanger.ac.uk, tgc@well.ox.ac.uk>, YY Teo <yy.teo@well.ox.ac.uk>
//
// $Id$
//

#define PI 3.141592654
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define log_mvt_t_density2_m 2

#define READ_CHAR_BUFFER (1024*1024*7)

#include <iostream.h>
#include <iomanip.h>
#include "./other_libraries/newmat11/newmatio.h"
#include "./other_libraries/newmat11/newmatap.h"

using namespace std;

extern "C" {
  int rmultinom_unif(int n);
  double runif();
  int rmultinom_unif(int n);
  void setall(long iseed1,long iseed2);
  float snorm(void);
  double pnorm(double x, double *mean, double *sd, int type);
  double dnorm(double x, double *mean, double *sd);
  double rnorm(double *mu, double *sd);
  double pt(double x, double df, int type);
  float genbet(float aa, float bb);
  double rgamma(double shape, double rate);
  float genchi(float df);
}

void matrix_exp(vector<vector<double> > *mat, double m);
void samp(int N, int M, int p, vector<int> *vec);
void samp1(int N, int M, int p, vector<int> *vec);
int countLines (istream& in);

int rmultinom(vector<double> *probs);
map<double, vector<int>, greater<double> >::iterator rmultinom(map<double, vector<int>, greater<double> > *m1, double total);

///////////////////////////////////////////////////////////////////////////
///// snp class which is used to format a snps data and for printing //////
///////////////////////////////////////////////////////////////////////////

class snp 
{
 public:
  double coord;
  string id;
  string rs;
  string allele1;
  string allele2;
  vector<ColumnVector> vec;
	
  // constructors/destructors	
  snp()
    {	
    }
  
  ~snp()
    {
    }
  
  void print() 
    {
      cout << rs << " " << id << endl;
      for(int j = 0; j < 10; j++) cout << "(" << vec[j](1) << ", " << vec[j](2) << ") ";
    cout << endl;
    }
  
};

/////////////////////////////////////////////////////////////////////////////////////////
////////// The data class contains all the variables, data and processing code //////////
/////////////////////////////////////////////////////////////////////////////////////////

class data
{
 public:
  int n_snp, n_ind, maxit, iter, cur_snp, niter, burnin, num_km, snp1, low, upp, best_start; 
  double thres, min_vec, max_vec;   
  vector<snp> dat;
  vector<double> strength;
  vector<double> contrast;
  vector<int> missing;
  vector<double> pert_score;
  bool print, calls, wga, probs, pert, not_in_pert;
  char *outfile; 
  
  // the calls 
  
  vector<int> ini_call; 
  vector<int> new_call; 
  vector<int> pert_call;
  vector<double> call_probs; // call probability matrix for a snp//
  vector<vector<int> > calls_all_snps; // calls for all snps //
  vector<vector<double> >  calls_all_probs; // call probs for all snps //
  
  // EM parameters
  
  vector<vector<double> > prs_mat; // initial parameter set of locations etc. // 
  vector<double> lambda;
  vector<vector<double> > sigma_mat;
  ColumnVector MU_AA, MU_AB, MU_BB, MU_NULL, dd;
  SymmetricMatrix covAA, covAB, covBB, covNULL ;	
  double muAA, muAB, muBB, muNULL, sdAA, sdAB, sdBB, sdNULL, maf;
  double mupert, sdpert;
  
  // data 
  
  ColumnVector tmpvec;
  vector<int> sex; // sex of individuals (1 = male, 0 = female)
  bool chrX; // flag for chrX SNPs
  vector<string> rs_id, illum_id, alleleA, alleleB;
  vector<int> bppos;
  vector<int> N;
  Matrix tmpmat;

  // storage
  
  SymmetricMatrix tmp;
  ColumnVector tmpcol;
  
  // constructors/destructors	
  data()
    {	
    }
  
  ~data()
    {
    }

//
// Split line into an array of words
//
void split_chars(char *s, vector <char*> &v) {
	char *p = strtok(s,"\t\n ");
	v.clear();
	while (p) {
		v.push_back(p);
		p = strtok(NULL,"\t\n ");
	}
}

  data(char *infile, char *xfile, int ll, int uu) {
		
		int total_number_of_snps = 0 ;
		char dummy[READ_CHAR_BUFFER] , *c , *d ;
		vector <char*> vc ;
		snp tmpsnp;
		tmpvec.ReSize(2);
		low = max ( ll , 1 ) ;
		upp = max ( -1 , uu ) ;
		if ( upp > 0 ) upp = max ( upp , low ) ;

		FILE *in_file = fopen ( infile , "r" ) ;
		if ( !in_file ) { cout << "Error : Could not open " << infile << endl ; exit(0) ; }
		
		// Header line
		fgets ( dummy , READ_CHAR_BUFFER , in_file ) ;
		split_chars ( dummy , vc ) ;
		n_ind = ( vc.size() - 3 ) / 2 ;
		cout << n_ind << " individuals" << endl ;
		
		// Read data
		while ( !feof ( in_file ) && ( upp == -1 || total_number_of_snps <= upp ) ) {
			fgets ( dummy , READ_CHAR_BUFFER , in_file ) ;

			split_chars ( dummy , vc ) ;
			if ( vc.size() < n_ind ) break ; // Last line
			
			total_number_of_snps++ ;
			if ( total_number_of_snps < low ) continue ;

			tmpsnp.rs = vc[0] ;
			tmpsnp.id = vc[1] ;
			tmpsnp.allele1 = vc[2] ;
			tmpsnp.vec.clear();
			for ( int a = 3 ; a < vc.size() ; a += 2 ) {
				tmpvec(1) = ( strcmp ( vc[a] , "NaN" ) == 0 ) ? -1.001 : atof ( vc[a] ) ;
				tmpvec(2) = ( strcmp ( vc[a+1] , "NaN" ) == 0 ) ? -1.001 : atof ( vc[a+1] ) ;
				tmpsnp.vec.push_back(tmpvec);
			}
			if ( tmpsnp.vec.size() != n_ind ) cout << tmpsnp.vec.size() ;
			dat.push_back(tmpsnp);
		}
		fclose ( in_file ) ;

		if ( upp == -1 ) upp = total_number_of_snps ;
		if ( upp > total_number_of_snps) {
			cout << "Warning: not all SNPs found; expecting at least " << upp << ", found " << total_number_of_snps << endl ;
			upp = total_number_of_snps ;
		}
		cout << "Total " << total_number_of_snps << " SNPs in file, " << low << " - " << upp << " processed." << endl ;
		n_snp = upp - low + 1;


		// Read chromosome X data, or not
		chrX = false;
		if ( xfile != NULL ) {
			FILE *x_file = fopen ( xfile , "r" ) ;
			while ( !feof ( x_file ) ) {
				fgets ( dummy , READ_CHAR_BUFFER , x_file ) ;
				if ( *dummy < 14 ) break ; // Last line
				sex.push_back ( atoi ( dummy ) ) ;
			}
			fclose ( x_file ) ;
			chrX = true;
		}

		// initialise some of the sizes of the data vectors 

		strength.resize(n_ind);
		contrast.resize(n_ind);
		ini_call.resize(n_ind);
		new_call.resize(n_ind);

		not_in_pert=true;

		vector<double> d1, d2;
		
		call_probs.resize((4*n_ind));
		d1.resize(3); for(int j = 0; j < 5; j++) prs_mat.push_back(d1);
		for(int j = 0; j < 5; j++) sigma_mat.push_back(d1); 

		lambda.resize(3);

		// crude initialisations - we fill in the other cells in the initialisation routine

		for(int j=0; j<5; j++) { for(int i = 0; i < 3; i++)  sigma_mat[j][i] = 0.1; }  

		prs_mat[0][0] = prs_mat[1][0] = prs_mat[2][0] = prs_mat[3][0] = -0.90000;
		prs_mat[0][2] = prs_mat[1][2]=  prs_mat[2][2] = prs_mat[4][2] =  0.90000;
		prs_mat[0][1] = 0.0; prs_mat[1][1]= -0.50000; prs_mat[2][1]= 0.5000;

		MU_AA.ReSize(2);   
		MU_AB.ReSize(2); 
		MU_BB.ReSize(2);
		MU_NULL.ReSize(2);

		covAA.ReSize(2);
		covAB.ReSize(2);
		covBB.ReSize(2);
		covNULL.ReSize(2);

		muNULL = 0.0;
		sdNULL = 100000.0;

		mupert = 0.0;
		sdpert = 0.05;
  }
  
  
  // class functions
  
  void illuminus() {

    if(pert) {
      pert_call.resize(n_ind);
      pert_score.resize(n_snp);
    }
    
    for(int i = 0; i < n_snp; i++) {
      cur_snp = i;
      if ( (i+1) % 100 == 0 ) cout << endl << "SNP " << i+1 << " " << flush;
      illum();
    }
    cout << endl;
    
    if(calls) output_illum_calls();
    if(probs) output_illum_probs();
    
    cout << "finished" << endl;
  }
  
  
  ////////////////////////////////
  /////// Output files ///////////
  ////////////////////////////////
  
  void output_illum_calls() {
    
    cout << "writing calls" << flush << endl;
    int i, j;
    
    ofstream ofile;
    string ss;
    
    ss = outfile;
    ss.append("_calls");
    
    ofile.open(ss.c_str());
    
    for(i = 0; i < n_snp; i++) {
      ofile << dat[i].rs << " " << dat[i].id << " " << dat[i].allele1 << " ";
      if(pert) ofile << setprecision(4) << pert_score[i] << " ";
      for(j = 0; j < n_ind; j++)  ofile << calls_all_snps[i][j] << " ";
      ofile << endl;
    }
    
    ofile.close();
  }
  
  void output_illum_probs() {
    
    cout << "writing probabilities" << flush << endl;
    int i, j;
    
    ofstream ofile;
    string ss;
    
    ss = outfile;
    ss.append("_probs");
    
    ofile.open(ss.c_str());
    
    for(i = 0; i < n_snp; i++) {
      ofile << dat[i].rs << " " << dat[i].id << " " << dat[i].allele1 << " ";   
      if(pert) ofile << setprecision(4) << pert_score[i] << " ";
      for(j = 0; j < (4*n_ind); j++) ofile << setprecision(4) << calls_all_probs[i][j] << " ";
      ofile << endl;
    }
    
    ofile.close();
    
  }
  
  /////////////////////////////////////////  
  /////////////////////////////////////////
  ///// THIS ROUTINE SETS IT ALL GOING ////
  /////////////////////////////////////////
  /////////////////////////////////////////
  
  
  void illum() {
    
    iter = 0;
    
    transform_intens(); // changes to contrast and strength;


    //// not all missing /////

    if( missing.size() < n_ind)
      {
	initialise_univ(); // sets up means and variances;
    
	choose_best_parameters(); // choosing the best starting position from the 5 possibles;
    
	initialise_mult();
    
	initialise_calls();
    
	while( update()  && iter < niter) {
	  iter++;
	}
    
	if(calls) calls_all_snps.push_back(ini_call); 
	if(probs) calls_all_probs.push_back(call_probs);	
	if(pert) perturbation(); 
      }

    //// all missing ////

    else
      {
	if(calls) 
	  {
	    for(int k = 0; k < n_ind; k++) ini_call[k] = 4;
	    calls_all_snps.push_back(ini_call);
	  }
	
	if(probs) 
	  {
	    for(int k=0; k < (4*n_ind); k++) 
	      {
		if( (k+1) % 4 == 0)  call_probs[k]=1.0;
		else call_probs[k]=0.0;
	      }
	    calls_all_probs.push_back(call_probs);
	  }

	if(pert)  pert_score.push_back(1.0);

      }
      
    // cout << iter << " " << flush;
    
  }
  
  
  ////////////////////////////////////
  // Transforming intensities ////////
  // /////////////////////////////////
  
  void transform_intens() {
    
    int j;
    double x1, y1;
    
    missing.clear();
    
    // cout << "transform intensities" << endl;
    for(j = 0; j < n_ind; j++) {
      x1 = (dat[cur_snp]).vec[j](1);
      y1 = (dat[cur_snp]).vec[j](2);
      if(  (x1>0.0 && y1 >=0.0) || (x1>=0.0 && y1>0.0) ) { 
        strength[j] = log( x1 + y1 );
        contrast[j] = (y1 - x1) / ( x1 + y1 );		
      }
      else {
        strength[j] = 0.0;
        contrast[j] = 0.0;
        missing.push_back(j);
      }
      // cout << x1 << " " << y1 << " " << strength[j] << " " << contrast[j] << endl;  
    } 
  }
  
  //////////////////////////////////////////////////////////////
  /////////////// initialise vectors and matrices //////////////
  ////////////////// for the univariate analysis ///////////////
  //////////////////////////////////////////////////////////////
  
  void initialise_univ() {
    
    // initialise prs_mat and sigma_mat//
   
    minmaxvec ( contrast , n_ind , min_vec , max_vec ) ;
  
    prs_mat[3][1] = 0.5 * ( max_vec + min_vec ); 
    prs_mat[4][1] = 0.5 * ( max_vec + min_vec ); 
    prs_mat[3][2] =  0.9 *  max_vec;
    prs_mat[4][0] =  0.9 *  min_vec;

    for(int j = 0; j < 2; j++)
      {
	for(int  i = 0; i < 3; i++)  { sigma_mat[j+3][i] = 0.05 * (max_vec - min_vec) ; }
      }
    
  }
  
  
  ///////////////////////////////////////////////////////////////
  ///////////// initialise for the multivariate setting /////////
  ///////////////////////////////////////////////////////////////
  
  void initialise_mult() {
    
    // mu // 
    
    MU_AA << muAA << 0.0;
    MU_AB << muAB << 0.0;
    MU_BB << muBB << 0.0;
    MU_NULL << 0.0 << 0.0;
    
    // covariance // 
    
    covAA.Row(1) << 0.05; covAA.Row(2) << 0.00 << 0.17;
    covAB.Row(1) << 0.05; covAB.Row(2) << 0.00 << 0.17;
    covBB.Row(1) << 0.05; covBB.Row(2) << 0.00 << 0.17;
    
    if(wga) { covNULL.Row(1) << 1000.0; covNULL.Row(2) << 0.00 << 1000.0; }
    else { covNULL.Row(1) << 100000.0; covNULL.Row(2) << 0.00 << 100000.0; }
    
    // lamda //
    
    for(int j = 0; j < 3; j++)  lambda[j]=1/3;
    
  }
  
  ///////////////////////////////////////////////////////////////
  /////////////// choose best starting parameters ///////////////
  ///////////////  from the five possible starts ////////////////
  ///////////////////////////////////////////////////////////////
  
  void choose_best_parameters() {
    
    vector<double> log_lik;
    vector<double>::iterator pos;
    int pos1;
    
    log_lik.resize(5);
    
    best_start=1;
    
    for(int j = 0; j < 5; j++) {
      muAA = prs_mat[j][0]; 
      muAB = prs_mat[j][1]; 
      muBB = prs_mat[j][2];
      sdAA = sigma_mat[j][0];
      sdAB = sigma_mat[j][1];
      sdBB = sigma_mat[j][2];
      
      initialise_calls();
      log_lik[j] = likelihood_ini();
      
    }
    
    pos = max_element( log_lik.begin() , log_lik.end() );		
    for(int k = 0; k < 5; k++) { if(log_lik[k] == *pos) pos1 = k; }
    
    // cout << "pos = " << pos1 << endl;
    muAA = prs_mat[pos1][0];
    muAB = prs_mat[pos1][1];
    muBB = prs_mat[pos1][2];
    sdAA = sigma_mat[pos1][0];
    sdAB = sigma_mat[pos1][1];
    sdBB = sigma_mat[pos1][2];
    
  }
  
  /////////////////////////////////////////////////////////////////////
  //////////////////// Initialise the calls ///////////////////////////
  /////////////////////////////////////////////////////////////////////
  
  void initialise_calls()  {
    
    vector<double> prob;
    vector<double>::iterator pos;
    double contr, dAA, dAB, dBB, pAA, pAB, pBB, sump, pNULL;
    int j, k, l, gc;
    
    double pnorm1 = pnorm( -1.0, &muAA, &sdAA, 1) ;
    double pnorm2 = pnorm( 1.0, &muAB, &sdAB, 0) - pnorm( -1.0, &muAB, &sdAB,0 ) ;
    double pnorm3 = pnorm(  1.0, &muBB, &sdBB, 0) ;
    
    prob.resize(4);
    for(j=0; j< n_ind; j++) {
      
      double contr = contrast[j];
      
      dAA = dnorm(contr, &muAA, &sdAA);
      dAB = dnorm(contr, &muAB, &sdAB);
      dBB = dnorm(contr, &muBB, &sdBB);
      pNULL = dnorm(contr, &muNULL, &sdNULL);
      
      pAA = dAA /  pnorm1; 
      pAB = dAB / pnorm2;
      pBB = dBB /  pnorm3;
      
      if(chrX && sex[j] == 1) prob[1] = 0.0; 
      if(chrX && sex[j] == 1) pAB = 0;

      sump = pAA + pAB + pBB + pNULL;
      
      prob[0] = pAA / sump;
      prob[1] = pAB / sump;
      prob[2] = pBB / sump;;
      prob[3] = pNULL / sump;
      
      pos = max_element( prob.begin() , prob.end() );
      gc = 4;
      for(k = 0; k < 4; k++) { if(prob[k] == *pos && prob[k] > thres) gc = k+1; }
      ini_call[j] = gc;
    }
    
    for(k = 0; k < missing.size(); k++) {
      l = missing[k]; 
      ini_call[l] = 4;
    }		
    
    new_call = ini_call; //.assign( ini_call.begin(), ini_call.end() );   
    
  }
  
  
  /////////////////////////////////////////////////////////////////////
  ///// Calculating the log likelihood for the univariate case ////////
  /////////////////////////////////////////////////////////////////////
  
  double likelihood_ini() {
    
    double loglike_AA, loglike_AB, loglike_BB, loglike_NULL;
    vector<double> cAA;
    vector<double> cAB;
    vector<double> cBB;
    vector<double> cNULL;
    double temp_mu, temp_sd, dn, pn, contr, total, pre_pnorm;
    int i, j;
    int size = 0;
    
    
    for(j = 0; j < n_ind; j++) {
      switch ( ini_call[j] ) {
      case 1 : cAA.push_back(contrast[j] ) ; break ;
      case 2 : cAB.push_back(  contrast[j] ); break ;
      case 3 : cBB.push_back(  contrast[j] ); break ;
      default : cNULL.push_back(contrast[j]); break ;
      }
    }
    
    // AA
    
    temp_mu = muAA; temp_sd = sdAA;
    size = cAA.size();
    loglike_AA = contr = dn = pn = total = 0.0;
    
    if( size >0 ) { 
      temp_mu = meanvec( cAA, size);
      if (size > 1 ) temp_sd = sqrt ( varvec( cAA, temp_mu, size) ) ;
      if ( temp_sd == 0.0 ) temp_sd = sdAA; 
    }
    
    pre_pnorm = pnorm(-1.0, &muAA, &temp_sd, 1 ) ;
    for( i = 0; i <  size; i++) {
      loglike_AA += log(dnorm(cAA[i], &temp_mu, &temp_sd)) / pre_pnorm;
    }
    
    // AB
    
    temp_mu = muAB; temp_sd = sdAB;
    size = cAB.size();
    loglike_AB = contr = dn = pn = 0.0;
    
    if( size >0 ) {
      temp_mu = meanvec( cAB, size);
      if (size > 1 ) temp_sd = sqrt ( varvec( cAB, temp_mu, size) ) ;
      if ( temp_sd == 0.0 ) temp_sd = sdAB;
    }
    
    pre_pnorm = pnorm(1.0, &muAB, &temp_sd, 0) - pnorm(-1.0, &muAB, &temp_sd, 0) ;
    for( i = 0; i <  size; i++) {
      dn = dnorm(cAB[i], &temp_mu, &temp_sd);
      loglike_AB += log(dn) / pre_pnorm;
    }
    
    // BB
    
    temp_mu = muBB; temp_sd = sdBB;
    size = cBB.size();
    loglike_BB = contr = dn = pn = 0.0;
    
    if( size >0 ) {
      temp_mu = meanvec( cBB, size);
      if (size > 1 ) temp_sd = sqrt ( varvec( cBB, temp_mu, size) ) ;
      if ( temp_sd == 0.0 ) temp_sd = sdBB;
    }
    
    pre_pnorm = pnorm(1.0, &muAA, &temp_sd, 0) ;
    for( i = 0; i <  size; i++) {
      dn = dnorm(cBB[i], &temp_mu, &temp_sd);
      loglike_BB += log(dn) / pre_pnorm;
    }
    
    // NULL
    
    temp_mu = muNULL; temp_sd = sdNULL;
    size = cNULL.size();
    loglike_NULL = contr = dn = pn = 0.0;
    
    if( size >0 ) {
      temp_mu = meanvec( cNULL, size);
      if (size > 1 ) temp_sd = sqrt ( varvec( cNULL, temp_mu, size) ) ;
      if ( temp_sd == 0.0 ) temp_sd = sdNULL;
    }
    
    for( i = 0; i <  size; i++) {
      loglike_NULL += log(dnorm(cNULL[i], &temp_mu, &temp_sd));
    }
    
    total = loglike_AA + loglike_AB + loglike_BB + loglike_NULL;
    return total;
  }
  
  /////////////////////////////////////////////////////////
  ///////////////////// Update routines ///////////////////
  /////////////////////////////////////////////////////////
  
  bool update() {
    
    maf = 0.0;
    
    int change=0; 
  
    // if(missing.size() == n_ind) return false;
 			  
    update_parameters();
    
    if(maf > 0.5) maf = 1 - maf;
    
    if( maf > 0.01 and ( covAA(1,1)  > covAB(1,1) or covBB(1,1)  > covAB(1,1))) change = 1;
    
    return update_prob(change);

  }
  
  ////////////////////////////
  // update probs and calls //
  ////////////////////////////
  
  
  bool update_prob(int change) {
   
    vector<double> prob;
    vector<double>::iterator pos;
    double contr, dAA, dAB, dBB, pAA, pAB, pBB, sump, pNULL, stdAA, stdAB1, stdAB2, stdBB;
    int j, k, l, gc, df;
    bool to_update;
    ColumnVector temp_mu;
    
    if(change == 1) df = 20;
    else df = 6;
    
    prob.resize(4);

    SymmetricMatrix covAB_i = covAB.i() ;
    SymmetricMatrix covAA_i = covAA.i() ;
    SymmetricMatrix covBB_i = covBB.i() ;
    SymmetricMatrix covNULL_i = covNULL.i() ;

    double lsd_AB = log(covAB.Determinant())/2 ;
    double lsd_AA = log(covAA.Determinant())/2 ;
    double lsd_BB = log(covBB.Determinant())/2 ;
    double lvd = -0.5* log( covNULL.Determinant() )  - log(2.0 * PI) ;

    stdAA = ( -1.0 - MU_AA(1) ) / sqrt( covAA(1,1) );
    stdAB1 = ( 1.0 - MU_AB(1) ) / sqrt( covAB(1,1) );
    stdAB2 = (-1.0 - MU_AB(1) ) / sqrt( covAB(1,1) );
    stdBB = (  1.0 - MU_BB(1) ) / sqrt( covBB(1,1) );
    
    double pt1 = pt( stdAA, df, 1) ;
    double pt2 = pt( stdAB1, 20, 0) - pt(stdAB2, 20, 0) ;
    double pt3 = pt( stdBB, df, 0) ;
    
    double lg1 = lgamma((log_mvt_t_density2_m + 20)/2)  - lgamma(20/2) ;
    double lg2 = lgamma((log_mvt_t_density2_m + df)/2)  - lgamma(df/2) ;
    
    for(j = 0; j < n_ind; j++) {
      
      temp_mu.ReSize(2);
      temp_mu << contrast[j] << strength[j];
      
      dAB = exp( log_mvt_t_density2( temp_mu - MU_AB, covAB_i, lg1-lsd_AB, 20));
      dAA = exp( log_mvt_t_density2( temp_mu - MU_AA, covAA_i, lg2-lsd_AA, df));
      dBB = exp( log_mvt_t_density2( temp_mu - MU_BB, covBB_i, lg2-lsd_BB, df));
      pNULL = exp( log_mvt_norm_density2( temp_mu - MU_NULL , covNULL_i, lvd));
  
      pAA = dAA / pt1;
      pAB = dAB / pt2;
      pBB = dBB / pt3 ;
      
      if( contrast[j] < MU_AB(1) ) pBB = 0.0;
      if( contrast[j] > MU_AB(1) ) pAA = 0.0;
      if( contrast[j] < MU_AA(1) or contrast[j] > MU_BB(1) ) pAB=0.0;
      
      pAA *= lambda[0];
      pAB *= lambda[1];
      pBB *= lambda[2];
      
      if(chrX && sex[j] == 1) prob[1] = 0.0; 
      if(chrX && sex[j] == 1) pAB = 0;
      sump = pAA + pAB + pBB + pNULL;
      
      prob[0] = pAA / sump;
      prob[1] = pAB / sump;
      prob[2] = pBB / sump;
      prob[3] = pNULL / sump;
      
      pos = max_element( prob.begin(), prob.end() );
      
      gc = 4; 
      for(k = 0; k < 4; k++) { if(prob[k] == *pos && prob[k] > thres) { gc = k+1; break;} }
      ini_call[j] = gc; 
      
      if(probs and not_in_pert) {
        call_probs[(j*4)] = prob[0];
        call_probs[(j*4+1)] = prob[1]; // heterozygous
        call_probs[(j*4+2)] = prob[2];	
        call_probs[(j*4+3)] = prob[3]; // null 
      }
      
    }
    
    for(k = 0; k < missing.size(); k++) {
      l = missing[k];
      ini_call[l] = 4; 
      if(probs)
	{
	  call_probs[(l*4)] = 0.0;
	  call_probs[(l*4+1)] = 0.0; // heterozygous
	  call_probs[(l*4+2)] = 0.0;	
	  call_probs[(l*4+3)] = 1.0; // null 
	}
    }    
    to_update = ini_call != new_call ;
    
    if(to_update) new_call = ini_call;//.assign( ini_call.begin(), ini_call.end() );
    
    return to_update;
  }
  
  
  //////////////////////////////////////////////////////
  ////////////// Update mean and covariance ////////////
  //////////////////////////////////////////////////////
  
  void update_parameters() {
    
    int j, l, nn, nn1;
    double s, nd;
    double pa, pb, totp;
    
    vector<double> ss1, ss2;
    ss1.resize(n_ind);
    ss2.resize(n_ind);
    
    maf = pa = pb = totp = 0.0;
    
    // AA
    for ( j = nn = 0 ; j < n_ind; j++) {
      if(ini_call[j] != 1) continue ;
      ss1[nn] = contrast[j] ;
      ss2[nn] = strength[j] ;
      nn++ ;
    }
   
    pa += 2*nn; 
    
    
    if(nn > 0) {               
      MU_AA(1) = meanvec(ss1, nn);               
      MU_AA(2) = meanvec(ss2, nn);               

    }                           

    if (nn > 4) {               
      covAA(1,1) = max( varvec(ss1, MU_AA(1), nn), 10e-6);
      covAA(2,2) = max( varvec(ss2, MU_AA(2), nn), 10e-6);               
      covAA(2,1) = covvec(ss1, ss2, MU_AA(1), MU_AA(2), nn);               
      
    }                        
    
    // BB                   
    
    for(nn = j = 0; j < n_ind; j++) {               
      if(ini_call[j] != 3) continue ;
      ss1[nn] = contrast[j] ;
      ss2[nn] = strength[j] ;
      nn++ ;
    }               
    pb += 2*nn; 
    
    if(nn > 0) {               
      
      MU_BB(1) = meanvec(ss1, nn);               
      MU_BB(2) = meanvec(ss2, nn);               
      
    }                         
    if (nn > 4) {               
      covBB(1,1) = max( varvec(ss1, MU_BB(1), nn), 10e-6) ;               
      covBB(2,2) = max( varvec(ss2, MU_BB(2), nn), 10e-6);               
      covBB(2,1) = covvec(ss1, ss2, MU_BB(1), MU_BB(2), nn);               
    }
    
    // AB                   
    
    for(nn = j = 0; j < n_ind; j++) {               
      if(ini_call[j] != 2) continue ;
      ss1[nn] = contrast[j] ;
      ss2[nn] = strength[j] ;
      nn++ ;
    }               
    
    pa +=nn;
    pb +=nn; 
    
    if(nn > 0) {               
      
      MU_AB(1) = meanvec(ss1, nn);               
      MU_AB(2) = meanvec(ss2, nn);                             
      
    }                         
    if (nn > 4) {               
      covAB(1,1) = max( varvec(ss1, MU_AB(1), nn), 10e-6);               
      covAB(2,2) = max( varvec(ss2, MU_AB(2), nn), 10e-6);               
      covAB(2,1) = covvec(ss1, ss2, MU_AB(1), MU_AB(2), nn);               
    }
    
    totp = pa + pb;

	if (chrX) {
		lambda[0] = 1.0;
		lambda[1] = 1.0;
		lambda[2] = 1.0;
	} else {
		lambda[0] = pa * pa / totp / totp;
		lambda[1] = 2 * pa * pb / totp / totp;
		lambda[2] = pb * pb / totp / totp;
	}
    
    maf = pa / totp; 
    
  }
  
  ///////////////////////////////////////////////
  //////////// Perturbation analysis ////////////
  ///////////////////////////////////////////////
  
  void perturbation() {
    
    int j=0;
    int l;
    double concord = 0.0;
    
    not_in_pert=false;
   
    for(int i=0; i<n_ind; i++) {

      strength[i] += rnorm(&mupert,&sdpert);
      contrast[i] += rnorm(&mupert,&sdpert);
    }
    
    for(int k = 0; k < missing.size(); k++) {
      l = missing[k];
      strength[l] = contrast[l] = 0.0; 
    }
    
    initialise_mult();
    
    pert_call = new_call = ini_call ;
         
    iter=0;
    while( update()  && iter < niter) { iter++; }
    
    for(int i=0; i<n_ind; i++) {
      if( pert_call[i] !=  ini_call[i] ) j++;
    }
    
    concord = (double) 1.0 - (double) j / (double) n_ind;
    pert_score[cur_snp] = concord; 
   
    not_in_pert=true;
    
  }
  
  //////////////////////////////////////////////
  //////// Summary statistics of vectors ///////
  //////////////////////////////////////////////
  
  inline void minmaxvec ( const vector<double> &vv , int n , double &min_vec , double &max_vec ) {
    min_vec = max_vec = vv[0] ;
    for ( register int a = 1 ; a < n ; a++ ) {
      if ( vv[a] < min_vec ) min_vec = vv[a] ;
      else if ( vv[a] > max_vec ) max_vec = vv[a] ;
    }
  }
  
  inline double meanvec(const vector<double> &vv, int n) {               
    
    double a = 0.0;  
    
    for(register int i = 0; i < n; i++) a += vv[i];               
    return a / n;               
    
  }               
  
  inline double varvec(const vector<double> &vv, double mu, int n) {               
    
    double a = 0.0;    

    for(register int i = 0; i < n; i++) a += (vv[i] - mu) * (vv[i] - mu);               
    return a / (n - 1);               
    
  }
  
  inline double covvec(const vector<double> &vv1, const vector<double> &vv2, double mu1, double mu2, int n) {               
    
    double a = 0.0;               
    for(register int i = 0; i < n; i++) a += (vv1[i] - mu1) * (vv2[i] - mu2);               
    return a / (n - 1);               
    
  }               
  
  /////////////////////////////////////
  ///////// log densities /////////////
  /////////////////////////////////////
  
  inline double log_mvt_norm_density2(const ColumnVector &xd, const SymmetricMatrix &V, double lvd) {               
    return lvd -0.5*(((xd.t()) * V) * xd).AsScalar();               
    
  }               
  
    
  inline double log_mvt_t_density2(const ColumnVector &xd, const SymmetricMatrix &sigma, double sig_det, int df)  {
    
    double a, distval;
    
    distval =   ( ( xd.t() * sigma ) * xd).AsScalar();
    a = sig_det
      - 0.5 * log_mvt_t_density2_m * log(PI * df)
      - 0.5 * (df + log_mvt_t_density2_m) * log(1 + distval/df);
    return a;
  }
  
  
  int countLines (istream& in)
    {
      return count(istreambuf_iterator<char>(in), istreambuf_iterator<char>(), '\n');
    }
  
  
};

