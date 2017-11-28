#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define Max(a,b) ((a)>(b)?(a):(b))

int it, itmax = 100;
int N, C, chunk_size;
int i, j, k, n, m, x, y;
double eps;
double maxeps = 0.1e-7;

double **A;

void relax();
void init();
double verify();

int main(int argc, char **argv) {
#ifndef _OPENMP
	fprintf(stderr, "Required OpenMP to run, check runtime parametres.\n");
	return 1;
#endif

	double time, sum;

	if (argc != 2)
	{
		fprintf(stderr, "Incorrect arguments number.\n");
		return 2;
	}

	sscanf(argv[1], "%d", &N);

	A = (double**)malloc(N * sizeof(double*));
	if (A == NULL)
	{
		fprintf(stderr, "Malloc error.\n");
		return 3;
	}
	for (i = 0; i < N; i++)
	{
		A[i] = (double*)malloc(N * sizeof(double));
		if (A[i] == NULL)
		{
			fprintf(stderr, "Malloc error.\n");
			return 3;
		}
	}
	init();

	printf("Matrix size: %d\n", N);
	time = omp_get_wtime();

	for (it = 1; it <= itmax; it++)
	{
		eps = 0.;
		relax();
		if (eps < maxeps) break;
	}

	time = omp_get_wtime() - time;
	printf("Time: %.20f\n", time);

	sum = verify();
	printf("Verify: %.20f\n", sum);

	for (i = 0; i < N; i++)
		free(A[i]);
	free(A);

	return 0;
}

void init()
{
	for (i = 0; i <= N - 1; i++)
		for (j = 0; j <= N - 1; j++)
		{
			if (i == 0 || i == N - 1 || j == 0 || j == N - 1)
				A[i][j] = 0.;
			else A[i][j] = (1. + i + j);
		}
}

void relax()
{
	for (i = 1; i <= N - 2; i++)
		for (j = 1; j <= N - 2; j++)
		{
			double e;
			e = A[i][j];
			A[i][j] = (A[i - 1][j] + A[i + 1][j] + A[i][j - 1] + A[i][j + 1]) / 4.;
			eps = Max(eps, fabs(e - A[i][j]));
		}
}

double verify()
{
	double s;

	s = 0.;
	for (i = 0; i <= N - 1; i++)
		for (j = 0; j <= N - 1; j++)
		{
			s = s + A[i][j] * (i + 1)*(j + 1) / (N * N);
		}
	return s;
}
