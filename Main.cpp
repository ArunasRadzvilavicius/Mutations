#include <iostream>
#include "Pop.h"

using namespace std;

int main(int argc, char* argv[]) {
	if (argc!=8) {cout << "Params N C R L K M mu" << endl; return 0;}
	int N = atoi(argv[1]);
	int C = atoi(argv[2]);
	double R = atof(argv[3]);
	double L = atof(argv[4]);
	int K = atoi(argv[5]);
	int M = atoi(argv[6]);
	double mu = atof(argv[7]);
	int B = 2;
	
	Population Pop(N, M, K, mu, C, R, L, B);
	Pop.Evolve(100000);
	//cout << Pop.NFixed() << endl;
	return 0;
}