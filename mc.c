// mc.c : Monte Carlo calculation on ising model with triangle lattice

#define DIM 2 // dimension
#define NN  6 // # of nearest-neighbors

#define L 2 // order of nearest-neighbor
#define N 7 // # of sites

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct Lattice {
	double s; // spin
} Lattice;

void Test(Lattice *lat) {
	FILE *f = fopen("test.txt", "w");
	int i, j, k, Ns = pow(2, N), Nt = 50;
	double e[Ns], m[Ns], sum, t, z, ee, mm;

	for(i=0; i<Ns; i++) {
		for(j=0; j<N; j++) lat[j].s = 2*((i >> j) & 0x01)-1;

		e[i] = m[i] = 0;
		for(j=0; j<N; j++) {
			sum = 0;
			for(k=0; k<N; k++) {
				if(k != j) sum += lat[k].s;
			}
			e[i] += -0.5 * lat[j].s * sum;
			m[i] += lat[j].s;
		}
	}

	for(i=1; i<Nt+1; i++) {
		t = 0.1 * i;	

		z = 0;
		for(j=0; j<Ns; j++)	z += exp(-e[j] / t);

		ee = mm = 0;
		for(j=0; j<Ns; j++) {
			ee += e[j] * exp(-e[j] / t) / z;
			mm += m[j] * exp(-e[j] / t) / z;
			//mm += fabs(m[j]) * exp(-e[j] / t) / z;
		}

		fprintf(f, "%f\t%f\t%f\n", t, ee, mm);
	}
	fclose(f);
}

int main(int argc, char *argv[]) {
	if(argc < 1) {
		printf("%s : Monte Carlo calculation on ising model with triangle lattice\n", argv[0]);
		exit(1);
	}

	Lattice lat[N];

	Test(lat);

	return 0;
}
