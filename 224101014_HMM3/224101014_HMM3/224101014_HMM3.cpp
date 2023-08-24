// 224101014_HMM3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#define N 5
#define M 32
#define T 85
#define MIN 1e-30

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

void backwardProcedure()
{
	long double sum = 0;
	for (int i = 1; i <= N; i++)
	{
		beta[T][i] = 1;
	}

	for (int t = T - 1; t >= 1; t--)
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				sum += (a[i][j] * b[j][o[t + 1]] * beta[t + 1][j]);
			}
			beta[t][i] = sum;
			sum = 0;
		}
	}
}

void calculateGamma()
{
	long double temp, sum;
	for (int t = 1; t <= T; t++)
	{
		sum = 0;
		for (int k = 1; k <= N; k++)
		{
			sum += (alpha[t][k] * beta[t][k]);
		}
		for (int i = 1; i <= N; i++)
		{
			temp = alpha[t][i] * beta[t][i];

			gamma[t][i] = (alpha[t][i] * beta[t][i]) / sum;
		}
	}
}

void viterbi()
{
	long double max, temp;
	int maxIndex;

	for (int i = 1; i <= N; i++)
	{
		delta[1][i] = pi[i] * b[i][o[1]];
		si[1][i] = 0;
	}

	for (int t = 2; t <= T; t++)
	{
		for (int j = 1; j <= N; j++)
		{
			max = 0;
			for (int i = 1; i <= N; i++)
			{
				temp = delta[t - 1][i] * a[i][j];
				if (max < temp)
				{
					max = temp;
					si[t][j] = i;
				}
			}
			delta[t][j] = max * b[j][o[t]];
		}
	}

	max = delta[T][1];
	maxIndex = 1;

	for (int i = 2; i <= N; i++)
	{
		if (max < delta[T][i])
		{
			max = delta[T][i];
			maxIndex = i;
		}
	}

	p_star = max;
	q_star[T] = maxIndex;

	for (int t = T - 1; t >= 1; t--)
	{
		q_star[t] = si[t + 1][q_star[t + 1]];
	}
}

void restimation()
{
	long double sum;

	for (int t = 1; t <= T - 1; t++)
	{
		sum = 0;
		for (int k = 1; k <= N; k++)
		{
			for (int m = 1; m <= N; m++)
			{
				sum += (alpha[t][k] * a[k][m] * b[m][o[t + 1]] * beta[t + 1][m]);
			}
		}
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				zeta[t][i][j] = ((alpha[t][i] * a[i][j] * b[j][o[t + 1]] * beta[t + 1][j]) / sum);
			}
		}
	}

	long double num, den;

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			num = den = 0;
			for (int t = 1; t <= T - 1; t++)
			{
				num += zeta[t][i][j];
				den += gamma[t][i];
			}
			a[i][j] = num / den;
		}
	}

	for (int j = 1; j <= N; j++)
	{
		for (int k = 1; k <= M; k++)
		{
			num = den = 0;
			for (int t = 1; t <= T; t++)
			{
				if (q_star[t] == j && o[t] == k)
				{
					num += gamma[t][j];
				}
				den += gamma[t][j];
			}
			b[j][k] = (num / den);
		}
	}
}

void reAdjustment(){
	long double max = 0, sum = 0;
	int index = 0;
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			sum += a[i][j];
			if(a[i][j] > max){
				max = a[i][j];
				index = j;
			}
		}
		if(sum <= 1){
			a[i][index] += (1 - sum);
		}
		else{
			a[i][index] -= (sum - 1);
		}
		sum = max = 0;
		index = 0;
	}
	
	index = 0;
	max = 0;
	sum =0;
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			if(b[i][j] < MIN){
				b[i][j] = MIN;
			}
			if(max < b[i][j]){
				max = b[i][j];
				index = j;
			}
			sum+=b[i][j];
		}
		if(sum > 1){
			b[i][index] -= (sum - 1);
		}
		else{
			b[i][index] += (1 - sum);
		}
		sum = index = max = 0;
	}
	
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

	fprintf(op, "\nBeta Matrix\n");
	fprintf(op, "================\n");
	for (int i = 1; i <= T; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(op, "%.32Le  ", beta[i][j]);
		}
		fprintf(op, "\n");
	}

	fprintf(op, "\nGamma Matrix\n");
	fprintf(op, "================\n");
	for (int i = 1; i <= T; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(op, "%.32Le  ", gamma[i][j]);
		}
		fprintf(op, "\n");
	}

	fprintf(op, "\nDelta Matrix\n");
	fprintf(op, "================\n");
	for (int i = 1; i <= T; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(op, "%.32Le  ", delta[i][j]);
		}
		fprintf(op, "\n");
	}

	fprintf(op, "\nSi Matrix\n");
	fprintf(op, "================\n");
	for (int i = 1; i <= T; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(op, "%d  ", si[i][j]);
		}
		fprintf(op, "\n");
	}

	fprintf(op, "\nQ star\n");
	fprintf(op, "================\n");
	for (int t = 1; t <= T; t++)
	{
		fprintf(op, "%d  ", q_star[t]);
	}

	fprintf(op, "\n");

	fprintf(op, "\nP star\n");
	fprintf(op, "================\n");
	fprintf(op, "%.32Le\n", p_star);

	fprintf(op, "\nZeta Matrix\n");
	fprintf(op, "================\n");

	for (int i = 1; i <= T; i++)
	{
		fprintf(op, "%d\n", i);
		for (int j = 1; j <= N; j++)
		{
			for (int k = 1; k < N; k++)
			{
				fprintf(op, "%.32Le", zeta[i][j][k]);
			}
			fprintf(op, "\n");
		}
		fprintf(op, "\n");
		fprintf(op, "\n");
	}

	fprintf(op, "\nA Matrix after adjustment\n");
	fprintf(op, "================\n");
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(op, "%.32Le  ", a[i][j]);
		}
		fprintf(op, "\n");
	}

	fprintf(op, "\nB Matrix after adjustment\n");
	fprintf(op, "================\n");
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			fprintf(op, "%.32Le  ", b[i][j]);
		}
		fprintf(op, "\n");
	}
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
		printf("oin not accesible");
	else
		printf("Reading oi values\n");

	for (int i = 1; i <= T; i++)
	{
		fscanf(oin, "%d", &o[i]);
	}

	op = fopen("output.txt", "w");

	if (!op)
		perror("output file not accesible");
	else
		printf("Writing output\n");

	forwardProcedure();
	backwardProcedure();
	calculateGamma();
	viterbi();
	restimation();
	reAdjustment();
	printOutput(op);
	fclose(op);
	getchar();
	return 0;
}

