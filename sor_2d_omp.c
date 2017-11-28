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
double **E;

void relax();
void init();
void init_eps();
double verify();

int main(int argc, char **argv) {
#ifndef _OPENMP
	fprintf(stderr, "Required OpenMP to run, check runtime parametres.\n");
	return 1;
#endif

	int max_threads, total_threads;
	double time, sum;
	max_threads = omp_get_max_threads();

	if (argc != 4)
	{
		fprintf(stderr, "Incorrect arguments number.\n");
		return 2;
	}

	sscanf(argv[1], "%d", &total_threads);
	sscanf(argv[2], "%d", &N);
	sscanf(argv[3], "%d", &chunk_size);
	C = (N - 2) / chunk_size;

	if (total_threads > max_threads)
	{
		fprintf(stderr, "Requested threads number (%d) is bigger than max thread for use (%d).\n", total_threads, max_threads);
		return 3;
	}

	A = (double**)malloc(N * sizeof(double*));
	if (A == NULL)
	{
		fprintf(stderr, "Malloc error.\n");
		return 4;
	}
	for (i = 0; i < N; i++)
	{
		A[i] = (double*)malloc(N * sizeof(double));
		if (A[i] == NULL)
		{
			fprintf(stderr, "Malloc error.\n");
			return 4;
		}
	}
	E = (double**)malloc(C * sizeof(double*));
	if (E == NULL)
	{
		fprintf(stderr, "Malloc error.\n");
		return 4;
	}
	for (i = 0; i < C; i++)
	{
		E[i] = (double*)malloc(C * sizeof(double));
		if (E[i] == NULL)
		{
			fprintf(stderr, "Malloc error.\n");
			return 4;
		}
	}

	printf("Matrix size: %d\n", N);
	omp_set_num_threads(total_threads);

	init();
	time = omp_get_wtime();
	for (it = 1; it <= itmax; it++)
	{
		eps = 0.;
		init_eps();
		relax();
		if (eps < maxeps) break;
	}
	time = omp_get_wtime() - time;
	printf("Time: %.20f\n", time);

	sum = verify();
	printf("Verify: %.20f\nThreads: %d\nChunk size: %d\nChunks amount: %d\n", sum, total_threads, chunk_size, C);

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

void init_eps()
{
	for (i = 0; i < C; i++)
		for (j = 0; j < C; j++)
			E[i][j] = 0.;
}

void relax()
{
	double e;
	for (k = 0; k < C; k++) {
#pragma omp parallel for shared(A, E, k) private(n, m, x, y, i, j, e)
		for (n = 0; n < k + 1; n++) {
			m = k - n;
			for (x = 0; x < chunk_size; x++)
				for (y = 0; y < chunk_size; y++)
				{
					i = x + 1 + n * chunk_size;
					j = y + 1 + m * chunk_size;
					e = A[i][j];
					A[i][j] = (A[i - 1][j] + A[i + 1][j] + A[i][j - 1] + A[i][j + 1]) / 4.;
					E[n][m] = Max(E[n][m], fabs(e - A[i][j]));
				}
		}

	}

	for (k = 0; k < C - 1; k++) {
#pragma omp parallel for shared(A, E, k) private(n, m, x, y, i, j, e)
		for (n = k + 1; n < C; n++) {
			m = C + k - n;
			for (x = 0; x < chunk_size; x++)
				for (y = 0; y < chunk_size; y++)
				{
					i = x + 1 + n * chunk_size;
					j = y + 1 + m * chunk_size;
					e = A[i][j];
					A[i][j] = (A[i - 1][j] + A[i + 1][j] + A[i][j - 1] + A[i][j + 1]) / 4.;
					E[n][m] = Max(E[n][m], fabs(e - A[i][j]));
				}
		}
	}

	for (i = 0; i < C; i++)
		for (j = 0; j < C; j++)
			if (eps < E[i][j]) eps = E[i][j];
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
