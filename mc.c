// mc.c : Monte Carlo calculation on ising model with triangle lattice

#define DIM 2    // dimension
#define NN  6    // # of nearest-neighbors

#define T_MAX 100 // maximum temperature
#define NT    128 // # of temperature interval

#define NORM(x) sqrt(x.c[0]*x.c[0] + x.c[1]*x.c[1])

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ran2.c"

typedef struct Vector {
	double c[DIM]; // coordinate
} Vector;

typedef struct Lattice {
	double s;               // spin
	double c[DIM];          // coordinate
	struct Lattice *nn[NN]; // addresses of nearest-neighbors
} Lattice;

long seed = -1;

void ReadVector(double *rad, Vector *r_lat, Vector *r_sup, Vector *r_sub) {
	FILE *f = fopen("input/vector.txt", "r");
	int i;

	fscanf(f, "%lf", rad);
	for(i=0; i<DIM;  i++) fscanf(f, "%lf%lf", &r_lat[i].c[0], &r_lat[i].c[1]);
	for(i=0; i<DIM;  i++) fscanf(f, "%lf%lf", &r_sup[i].c[0], &r_sup[i].c[1]);
	for(i=0; i<NN+1; i++) fscanf(f, "%lf%lf", &r_sub[i].c[0], &r_sub[i].c[1]);

	fclose(f);
}

void InitLat(Lattice *lat, double rad, Vector *r_sup, Vector *r_sub) {
	int i, j, k, cnt;
	int n0, n1, min = -4, max = 4;
	Vector r;

	for(i=0; i<NN+1; i++) {
		lat[i].s = ran2(&seed) > 0.5 ? 1 : -1;

		for(j=0; j<DIM; j++) lat[i].c[j] = r_sub[i].c[j];

		cnt = 0;
		for(n0=min; n0<max; n0++) {
			for(n1=min; n1<max; n1++) {
				for(j=0; j<NN+1; j++) {
					for(k=0; k<DIM; k++) r.c[k] = n0 * r_sup[0].c[k] + n1 * r_sup[1].c[k] + r_sub[j].c[k] - r_sub[i].c[k];
					if(fabs(NORM(r) - rad) < 1e-6) {
						lat[i].nn[cnt] = &lat[j];
						cnt++;
					}
				}
			}
		}

		if(cnt > NN) {
			printf("%d neighbors are expected but %d obtained", NN, cnt);
			exit(1);
		}
	}
}

double CalcEnergy(Lattice *lat) {
	int i, j;
	double e = 0;

	for(i=0; i<NN+1; i++) {
		for(j=0; j<NN; j++) {
			e += lat[i].s * lat[i].nn[j]->s;
		}
	}

	return -0.5 * e;
}

double CalcMagnetization(Lattice *lat) {
	int i;
	double m = 0;

	for(i=0; i<NN+1; i++) {
		m += lat[i].s;
	}

	return m;
}

void Test(Lattice *lat, double *T) {
	FILE *f;
	char fn[64];

	sprintf(fn, "output/test.txt");
	f = fopen(fn, "w");

	int i, j, Ns = pow(2, NN+1);
	double e[Ns], m[Ns], Z, E, M;

	time_t t0 = time(NULL);

	for(i=0; i<Ns; i++) {
		for(j=0; j<NN+1; j++) lat[j].s = 2*((i >> j) & 0x01)-1;

		e[i] = CalcEnergy(lat);
		m[i] = CalcMagnetization(lat);
	}

	for(i=0; i<NT; i++) {
		Z = E = M = 0;
		for(j=0; j<Ns; j++) Z += exp(-e[j] / T[i]);
		for(j=0; j<Ns; j++) {
			E +=       e[j] * exp(-e[j] / T[i]) / Z;
			M += fabs(m[j]) * exp(-e[j] / T[i]) / Z;
		}

		fprintf(f, "%f\t%f\t%f\n", T[i], E, M);
	}
	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
}

void MonteCarlo(Lattice *lat, double *T, int Nmc0) {
	FILE *f;
	char fn[64];

	sprintf(fn, "output/mc%d.txt", Nmc0);
	f = fopen(fn, "w");

	int i, j, k, n, Nmc = pow(10, Nmc0);
	double E0, E1, E, M;

	time_t t0 = time(NULL);

	for(i=0; i<NT; i++) {
		E = M = 0;
		for(j=0; j<Nmc; j++) {
			for(k=0; k<NN+1; k++) {
				n = (int)((NN+1) * ran2(&seed)) % (NN+1);			

				E0 = CalcEnergy(lat);
				lat[n].s *= -1;
				E1 = CalcEnergy(lat);

				if(E1 - E0 > 0 && ran2(&seed) > exp((E0 - E1) / T[i])) lat[n].s *= -1;
			}
			E += CalcEnergy(lat);
			M += fabs(CalcMagnetization(lat));
		}
		fprintf(f, "%f\t%f\t%f\n", T[i], E/Nmc, M/Nmc);
	}
	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
}

int main(int argc, char *argv[]) {
	if(argc < 1) {
		printf("%s <test/mc> <(mc)Nmc> : Monte Carlo calculation on ising model with triangle lattice\n", argv[0]);
		exit(1);
	}

	double rad;
	Vector r_lat[DIM], r_sup[DIM], r_sub[NN+1];
	ReadVector(&rad, r_lat, r_sup, r_sub);

	Lattice lat[NN+1];
	InitLat(lat, rad, r_sup, r_sub);

	int i;
	double T[NT];
	for(i=0; i<NT; i++) T[i] = T_MAX - T_MAX * ((double)i / NT);

	if(strstr(argv[1], "t")) Test(lat, T);
	else                     MonteCarlo(lat, T, atoi(argv[2]));

	return 0;
}
