#include <iostream>
#include "Pop.h"

using namespace std;

int main(int argc, char* argv[]) {
	int N = atoi(argv[1]);
	int C = atoi(argv[2]);
	double R = atof(argv[3]);
	double L = atof(argv[4]);
	Population Pop(N, 20, 100, 1.0, C, R, L);
	Pop.Evolve(25000);
	return 0;
}