#include <iostream>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <numeric>
#include <cmath>
#include <assert.h>
#include "Pop.h"

using namespace std;

Population::Population(){
}

Population::Population(int Nv, int Mv, int Kv, double muv, int Cv, double Rv, double Lv) {
	    
	N = Nv;
	M = Mv;
	K = Kv;
	C = Cv;				// size of a segregtaing unit in mtDNA's
	R = Rv;				// mixing rate at the lowest level - number of alleles randomly exchanged
	L = Lv;				// mixing rate at the higher level - number of llu randomly exchanged
	assert(M%C==0);
	assert(C<=M);
	assert(R<K);
	double mu0 = muv;	// mutation rate per genome per generation
	mu = mu0 / (1);	// mutation rate per site per generation
	s = 0.05;			// fitness contribution at the lowest level
	c = 1.0;			// fitness contrbution one level up
	Pop = vector<cell> (N);
	cell cl(M);
	group DNA(K);
	fill(DNA.begin(), DNA.end(), 0);
	fill(cl.begin(), cl.end(), DNA);
	fill(Pop.begin(), Pop.end(), cl);
	
	const gsl_rng_type *T;
	T = gsl_rng_mt19937;
	rng = gsl_rng_alloc(T);
	gsl_rng_set (rng, time(NULL));	
	
}

void Population::Evolve(int Gv){
	for (int i=0;i<Gv;i++) {
		//cout <<MeanMut()<<" "<<BestClass()<<" "<<NFixed()<<" "<<HLLC()<<" "<<HLLC1()<<endl;
		Mutate();
		Mix();
		Recombine();
		Division(C);
		Selection();
	}
}

void Population::Mutate(){
	for (int i=0;i<N;i++){
		int dk = gsl_ran_poisson(rng, mu);
		for (int k=0;k<dk;k++){
			Pop[i][gsl_rng_uniform_int(rng,M)][gsl_rng_uniform_int(rng,K)]+=1;
		}
	}
}

void Population::Selection(){
	vector<double> Fits(N);
	for (int i=0;i<N;i++){
		double mfit = 0.0;
		for (int j=0;j<M;j++){
			int nmut = accumulate(Pop[i][j].begin(), Pop[i][j].end(), 0);
			mfit = mfit + pow(1-s, nmut);
		}
		Fits[i] = 1.0 - c*pow(1.0 - mfit/M,1.0);
	}
	vector<int> S = WS(Fits, N);
	vector<cell> nPop(0);
	for (int i=0;i<N;i++){
		nPop.push_back(Pop[S[i]]);
	}
	Pop = vector<cell>(0);
	copy(nPop.begin(), nPop.end(), back_inserter(Pop));
	assert(Pop.size()==N);
}

void Population::Recombine(){
	for (int i=0;i<N;i++){
		int dr = gsl_ran_poisson(rng, R);
		for (int r=0;r<dr;r++){
			int k = gsl_rng_uniform_int(rng,K);
			int m1 = gsl_rng_uniform_int(rng,M);
			int m2 = gsl_rng_uniform_int(rng,M);
			Pop[i][m2][k] = Pop[i][m1][k];
		}		
	}
}

void Population::Mix(){
	cell cl1;
	cell cl2;
	cell sc1;
	cell sc2;
	random_shuffle(Pop.begin(), Pop.end());
	int i = 0;
	if (L>0){
		while (i<N){
			int l = gsl_ran_poisson (rng, L);
			if (l>M) {l=M;}
			cl1 = Pop[i];
			cl2 = Pop[i+1];
			random_shuffle(cl1.begin(), cl1.end());
			random_shuffle(cl2.begin(), cl2.end());
			sc1 = cell(cl1.begin(), cl1.begin()+l);
			sc2 = cell(cl2.begin(), cl2.begin()+l);
			copy(sc1.begin(), sc1.end(), cl2.begin());
			copy(sc2.begin(), sc2.end(), cl1.begin());
			Pop[i] = cl1;
			Pop[i+1] = cl2;
			i = i + 2;
		}
	}
}

void Population::Division(int C){
	vector<cell> nPop(N);
	int CN = M/C;					// number of clusters of size C
	cell cl(M);
	cell ncl;
	cell cluster(C);
	clustergroup Ccl(0);
	clustergroup DCcl(0);			// dublicated cluster cell
	cell Dcl(0);					// duplicated cell, unclustered
	for (int i=0;i<N;i++){
		Dcl = cell(0); 
		DCcl = clustergroup(0);
		cl = Pop[i];
		random_shuffle(cl.begin(), cl.end());
		Ccl = clustergroup (0);
		for (int j=0;j<CN;j++){		// make clusters within the clustercell
			cluster = cell( cl.begin()+j*C, cl.begin()+(j+1)*C );
			assert(cluster.size()==C);
			Ccl.push_back(cluster);
		}
		// duplicate
		copy(Ccl.begin(), Ccl.end(), back_inserter(DCcl));
		copy(Ccl.begin(), Ccl.end(), back_inserter(DCcl));
		random_shuffle(DCcl.begin(), DCcl.end());
		// divide
		Ccl = clustergroup(0);
		copy(DCcl.begin(), DCcl.begin()+CN, back_inserter(Ccl));
		//	uncluster
		ncl = cell(0);
		for (int k=0;k<CN;k++){
			cluster = Ccl[k];
			assert(cluster.size()==C);
			copy(Ccl[k].begin(), Ccl[k].end(), back_inserter(ncl));
		}		
		nPop[i] = ncl;
	}
	Pop = nPop;
}

vector<int> Population::WS(vector<double> weights, int Rw){
	double total = accumulate(weights.begin(),weights.end(), 0.0);
	vector<int> Lw (0);
	int i = 0;
	double w = weights[0];
	while (Rw>0){
		double r = gsl_rng_uniform (rng);
		double x = total*(1.0 - pow(r, 1.0/Rw));
		total = total - x;
		while (x>w){
			x = x-w;
			i += 1;
			w = weights[i];
		}
		w -= x;
		Rw -= 1 ;
		Lw.push_back(i);
	}
	return Lw;
}

void Population::PV(vector<int> v){
	for (int i=0;i<v.size();i++){
		cout << v[i] << " ";
	}
	cout << endl;
}

int Population::BestClass(){
	vector<int> Bests(0);
	for (int i=0;i<N;i++){
		for (int j=0;j<M;j++){
			Bests.push_back(accumulate(Pop[i][j].begin(), Pop[i][j].end(), 0));
		}
	}
	return *min_element(Bests.begin(), Bests.end());
}

int Population::NFixed(){
	vector<int> minsites(K,0);
	for (int k=0;k<K;k++){
		vector<int> mutsites(0);
		for (int i=0;i<N;i++){
			for (int j=0;j<M;j++){
				mutsites.push_back(Pop[i][j][k]);
			}
		}
		minsites[k] = *min_element(mutsites.begin(), mutsites.end());
	}
	return accumulate(minsites.begin(), minsites.end(), 0);
}

double Population::HLLC(){
	int mutllc = BestClass();
	vector<group> LLC(0);
	// isolate the LLC
	for (int i=0;i<N;i++){
		for (int j=0;j<M;j++){
			if (accumulate(Pop[i][j].begin(),Pop[i][j].end(),0)==mutllc){
				LLC.push_back(Pop[i][j]);
			}
		}
	}
	int Nllc = LLC.size();
	double H = 0;
	for (int k=0;k<K;k++){
		double sumx = 0;
		vector<int> mutsk(0);
		for (int i=0;i<Nllc;i++){
			mutsk.push_back(LLC[i][k]);
		}
		for (int x=0;x<=mutllc;x++){
			int fx = count(mutsk.begin(),mutsk.end(),x);
			double f = fx/(double)Nllc;
			sumx += pow(f,2);
		}
		H += 1.0 - sumx;
	}
	return H/K;
}

double Population::HLLC1(){
	int mutllc = BestClass()+1;
	vector<group> LLC(0);
	// isolate the LLC
	for (int i=0;i<N;i++){
		for (int j=0;j<M;j++){
			if (accumulate(Pop[i][j].begin(),Pop[i][j].end(),0)==mutllc){
				LLC.push_back(Pop[i][j]);
			}
		}
	}
	int Nllc = LLC.size();
	double H = 0;
	for (int k=0;k<K;k++){
		double sumx = 0;
		vector<int> mutsk(0);
		for (int i=0;i<Nllc;i++){
			mutsk.push_back(LLC[i][k]);
		}
		for (int x=0;x<=mutllc;x++){
			int fx = count(mutsk.begin(),mutsk.end(),x);
			double f = fx/(double)Nllc;
			sumx += pow(f,2);
		}
		H += 1.0 - sumx;
	}
	return H/K;
}

double Population::MeanMut(){
	vector<int> MMut(0);
	for (int i=0;i<N;i++){
		for (int j=0;j<M;j++){
			MMut.push_back(accumulate(Pop[i][j].begin(), Pop[i][j].end(), 0));
		}
	}
	return accumulate(MMut.begin(),MMut.end(),0.0)/(N*M);
}