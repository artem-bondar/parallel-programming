#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define Max(a,b) ((a)>(b)?(a):(b))

int it, itmax = 100;
int N, mpi_size, chunk_size, rank;
int i, j, k, l, n, m, x, y;
double eps, local_eps;
double maxeps = 0.1e-7;
double sum = 0, local_sum = 0;

double **A;
double *msg;

void relax();
void init();
void verify();

int main(int argc, char **argv) {
	int error;
	double time;

	if ((error = MPI_Init(&argc, &argv))) {
		fprintf(stderr, "MPI startup error!\n");
		MPI_Abort(MPI_COMM_WORLD, error);
		return 1;
	}

	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc != 2)
	{
		fprintf(stderr, "Incorrect arguments number.\n");
		MPI_Abort(MPI_COMM_WORLD, 2);
		return 2;
	}

	sscanf(argv[1], "%d", &N);
	chunk_size = (N - 2) / mpi_size;

	A = (double**)malloc(N * sizeof(double*));
	if (A == NULL)
	{
		fprintf(stderr, "Malloc error.\n");
		MPI_Abort(MPI_COMM_WORLD, 3);
		return 3;
	}
	for (i = 0; i < N; i++)
	{
		A[i] = (double*)malloc(N * sizeof(double));
		if (A[i] == NULL)
		{
			fprintf(stderr, "Malloc error.\n");
			MPI_Abort(MPI_COMM_WORLD, 3);
			return 3;
		}
	}

	msg = (double*)malloc(chunk_size * sizeof(double));
	if (msg == NULL)
	{
		fprintf(stderr, "Malloc error.\n");
		MPI_Abort(MPI_COMM_WORLD, 3);
		return 3;
	}

	init();
	MPI_Barrier(MPI_COMM_WORLD);
	if (!rank)
	{
		time = MPI_Wtime();
		printf("Matrix size: %d\n", N);
	}
	for (it = 1; it <= itmax; it++)
	{
		local_eps = eps = 0.;
		relax();
		MPI_Allreduce(&local_eps, &eps, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if (eps < maxeps) break;
	}
	if (!rank)
	{
		time = MPI_Wtime() - time;
		printf("Time: %.20f\n", time);
	}

	verify();
	MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (!rank)
		printf("Verify: %.20f\nProcesses: %d\nChunk size: %d\n", sum, mpi_size, chunk_size);

	for (i = 0; i < N; i++)
		free(A[i]);
	free(A);

	MPI_Finalize();
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
	double e;
	for (k = 0; k < mpi_size; k++) {
		for (n = 0; n < k + 1; n++) {
			m = k - n;
			if (rank == m)
			{
				if (m != 0)
				{
					MPI_Recv(msg, chunk_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					j = m * chunk_size;
					for (l = 0; l < chunk_size; l++)
					{
						i = l + 1 + n * chunk_size;
						A[i][j] = msg[l];
					}
				}
				for (x = 0; x < chunk_size; x++)
					for (y = 0; y < chunk_size; y++)
					{
						i = x + 1 + n * chunk_size;
						j = y + 1 + m * chunk_size;
						e = A[i][j];
						A[i][j] = (A[i - 1][j] + A[i + 1][j] + A[i][j - 1] + A[i][j + 1]) / 4.;
						local_eps = Max(local_eps, fabs(e - A[i][j]));
					}
				if (k != mpi_size - 1)
				{
					if (m != mpi_size - 1)
					{
						j = (m + 1) * chunk_size;
						for (l = 0; l < chunk_size; l++)
						{
							i = l + 1 + n * chunk_size;
							msg[l] = A[i][j];
						}
						MPI_Send(msg, chunk_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
					}
				}
				else
				{
					if (m != 0)
					{
						i = (n + 1) * chunk_size;
						for (l = 0; l < chunk_size; l++)
						{
							j = l + 1 + m * chunk_size;
							msg[l] = A[i][j];
						}
						MPI_Send(msg, chunk_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
					}
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	for (k = 0; k < mpi_size - 1; k++) {
		for (n = k + 1; n < mpi_size; n++) {
			m = mpi_size + k - n;
			if (rank == mpi_size - 1 - n)
			{
				MPI_Recv(msg, chunk_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				i = n * chunk_size;
				for (l = 0; l < chunk_size; l++)
				{
					j = l + 1 + m * chunk_size;
					A[i][j] = msg[l];
				}
				for (x = 0; x < chunk_size; x++)
					for (y = 0; y < chunk_size; y++)
					{
						i = x + 1 + n * chunk_size;
						j = y + 1 + m * chunk_size;
						e = A[i][j];
						A[i][j] = (A[i - 1][j] + A[i + 1][j] + A[i][j - 1] + A[i][j + 1]) / 4.;
						local_eps = Max(local_eps, fabs(e - A[i][j]));
					}
				if (n != mpi_size - 1)
				{
					i = (n + 1) * chunk_size;
					for (l = 0; l < chunk_size; l++)
					{
						j = l + 1 + m * chunk_size;
						msg[l] = A[i][j];
					}
					MPI_Send(msg, chunk_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	for (k = 0; k < mpi_size; k++) {
		for (n = 0; n < k + 1; n++) {
			m = k - n;
			if (rank == m)
			{
				if (m != 0)
				{
					j = m * chunk_size + 1;
					for (l = 0; l < chunk_size; l++)
					{
						i = l + 1 + n * chunk_size;
						msg[l] = A[i][j];
					}
					MPI_Send(msg, chunk_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
				}
				if (k != mpi_size - 1)
				{
					if (m != mpi_size - 1)
					{
						MPI_Recv(msg, chunk_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						j = (m + 1) * chunk_size + 1;
						for (l = 0; l < chunk_size; l++)
						{
							i = l + 1 + n * chunk_size;
							A[i][j] = msg[l];
						}
					}
				}
				else
				{
					if (m != 0)
					{
						MPI_Recv(msg, chunk_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						i = (n + 1) * chunk_size + 1;
						for (l = 0; l < chunk_size; l++)
						{
							j = l + 1 + m * chunk_size;
							A[i][j] = msg[l];
						}
					}
				}
			}
		}
	}
	for (k = 0; k < mpi_size - 1; k++) {
		for (n = k + 1; n < mpi_size; n++) {
			m = mpi_size + k - n;
			if (rank == mpi_size - 1 - n)
			{
				i = n * chunk_size + 1;
				for (l = 0; l < chunk_size; l++)
				{
					j = l + 1 + m * chunk_size;
					msg[l] = A[i][j];
				}
				MPI_Send(msg, chunk_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
				if (n != mpi_size - 1)
				{
					MPI_Recv(msg, chunk_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					i = (n + 1) * chunk_size + 1;
					for (l = 0; l < chunk_size; l++)
					{
						j = l + 1 + m * chunk_size;
						A[i][j] = msg[l];
					}
				}
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void verify()
{
	for (k = 0; k < mpi_size; k++) {
		for (n = 0; n < k + 1; n++) {
			m = k - n;
			if (rank == m)
			{
				for (x = 0; x < chunk_size; x++)
					for (y = 0; y < chunk_size; y++)
					{
						i = x + 1 + n * chunk_size;
						j = y + 1 + m * chunk_size;
						local_sum += A[i][j] * (i + 1)*(j + 1) / (N * N);
					}
			}
		}
	}
	for (k = 0; k < mpi_size - 1; k++) {
		for (n = k + 1; n < mpi_size; n++) {
			m = mpi_size + k - n;
			if (rank == mpi_size - 1 - n)
			{
				for (x = 0; x < chunk_size; x++)
					for (y = 0; y < chunk_size; y++)
					{
						i = x + 1 + n * chunk_size;
						j = y + 1 + m * chunk_size;
						local_sum = local_sum + A[i][j] * (i + 1)*(j + 1) / (N * N);
					}
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
