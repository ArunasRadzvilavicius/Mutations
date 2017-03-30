#include <iostream>
#include "Pop.h"

using namespace std;

int main(int argc, char* argv[]) {
	int N = atoi(argv[1]);
	int C = atoi(argv[2]);
	double R = atof(argv[3]);
	double L = atof(argv[4]);
	int K = atoi(argv[5]);
	Population Pop(N, 20, K, 0.1, C, R, L);
	Pop.Evolve(100000);
	//cout << Pop.NFixed() << endl;
	return 0;
}