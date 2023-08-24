// 224101014_HMM1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdafx.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#define N 5
#define M 32
#define T 85

long double a[N + 1][N + 1], b[N + 1][M + 1], pi[N + 1];
long double alpha[T + 1][N + 1];
long double beta[T + 1][N + 1];
long double gamma[T + 1][N + 1];
long double delta[T + 1][N + 1];
int si[T + 1][N + 1];
long double theta[T + 1][N + 1];
long double zeta[T + 1][N + 1][N + 1];
long double p_star;
int q_star[T + 1];
long double p;
int o[T + 1];

void initialisation()
{
	for (int i = 0; i <= T; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			alpha[i][j] = 0;
			beta[i][j] = 0;
			gamma[i][j] = 0;
			theta[i][j] = 0;
			si[i][j] = 0;
		}
	}

	for (int i = 0; i <= T; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int k = 0; k <= N; k++)
			{
				zeta[i][j][k] = 0;
			}
		}
	}

	for (int i = 0; i <= T; i++)
	{
		q_star[i] = 0;
	}
}

void forwardProcedure()
{
	long double sum = 0;

	for (int i = 1; i <= N; i++)
	{
		alpha[1][i] = pi[i] * b[i][o[1]];
	}

	for (int t = 1; t <= T - 1; t++)
	{
		for (int j = 1; j <= N; j++)
		{
			for (int i = 1; i <= N; i++)
			{
				sum += (alpha[t][i] * a[i][j]);
			}
			alpha[t + 1][j] = sum * b[j][o[t + 1]];
			sum = 0;
		}
	}

	sum = 0;

	for (int m = 1; m <= N; m++)
	{
		sum += alpha[T][m];
	}
	p = sum;
}

void printOutput(FILE *op)
{
	for (int i = 1; i <= T; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(op, "%.32Le  ", alpha[i][j]);
		}
		fprintf(op, "\n");
	}

	fprintf(op, "\np(O|lamda) = %.32Le\n", p);
}

int _tmain(int argc, _TCHAR* argv[])
{	
	FILE *ain, *bin, *pin, *oin, *op;
	ain = fopen("a.txt", "r");

	if (!ain)
		perror("a not accesible");
	else
		printf("Reading a values\n");

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fscanf(ain, "%Le", &a[i][j]);
		}
	}

	bin = fopen("b.txt", "r");

	if (!bin)
		perror("b not accesible");
	else
		printf("Reading b values\n");

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			fscanf(bin, "%Le", &b[i][j]);
		}
	}

	pin = fopen("pi.txt", "r");

	if (!pin)
		perror("pi not accesible");
	else
		printf("Reading pi values\n");

	for (int i = 1; i <= N; i++)
	{
		fscanf(pin, "%Le", &pi[i]);
	}

	oin = fopen("o.txt", "r");

	if (!oin)
		perror("oin not accesible");
	else
		printf("Reading oi values\n");

	for (int i = 1; i <= T; i++)
	{
		fscanf(oin, "%d", &o[i]);
	}

	for (int i = 1; i <= N; i++)
	{
		fscanf(pin, "%Le", &pi[i]);
	}

	op = fopen("output.txt", "w");

	if (!op)
		perror("output file not accesible");
	else
		printf("Writing output\n");

	forwardProcedure();
	printOutput(op);

	fclose(op);
	return 0;
}

