#ifndef POP_H_INCLUDED
#define POP_H_INCLUDED
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
using namespace std;

typedef vector<int> group;
typedef vector<group> cell;
typedef vector<cell> clustergroup;
class Population {
	public:
		Population();
		Population(int,int,int,double,int,double,double);
		void Evolve(int);
		void Mutate();
		void Selection();
		void Division(int);
		void Mix();
		void Recombine();
		vector<int> WS(vector<double>, int);
		vector<int> WSS(vector<double>, int);
		void PV(vector<int>);
		int BestClass();
		int NFixed();
		double HLLC();
		double HLLC1();
		double MeanMut();
		
		gsl_rng *rng;
		vector<cell> Pop;
		int N; // population size
		int M; // cell size
		int K; // number of mito loci
		int C; // segregating unit
		double R; // recombination rate
		double L; // mixing rate 2
		double mu;
		double s;
		double c;
};

#endif 
