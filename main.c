#include <stdio.h>
#include <stdlib.h>
#include "TSpack.h"

/** Main entry point. */
int main(int argc, char **argv)
{                                           /* driver program for the TSPACK code */
	int    i, j, value, correct, uheur;
	int    d, n, lb, ub0, ub;
	float  timeLimit;
	char   fname[40];
	FILE   *file;

	int    **w, **x;
	int    *W, *b;

	if (argc < 3) {
	   printf("\n Usage: %s [totalTime] [InputFile]\n",argv[0]);
	   exit(0);
	}

	sscanf(argv[1], "%f", &timeLimit);
	sscanf(argv[2], "%s", fname);
	if (argc == 4) {
       uheur = sscanf(argv[3], "%d", &uheur);
    } else {
       uheur = 1;
    }

	file = fopen(fname, "r");
	fscanf(file, "%d", &d);
    fscanf(file, "%d", &n);

	/* memory allocation */
	W = (int*)calloc(d, sizeof(int));
	b = (int*)calloc(n, sizeof(int));
	
	w = (int**)calloc(d, sizeof(int*));
	x = (int**)calloc(d, sizeof(int*));
	for (i = 0; i < d; i++) {
	    w[i] = (int*)calloc(n, sizeof(int));
	    x[i] = (int*)calloc(n, sizeof(int));
	}

	/* get the current instance */
	for (i = 0; i < d; i++) {
	    fscanf(file, "%d", &value);
	    W[i] = value;
	}
	for (i = 0; i < n; i++) {
	    for (j = 0; j < d; j++) {
	        fscanf(file, "%d", &value);
	        w[j][i] = value;
	    }
	}
	fclose(file);

	/* compute an initial lower bound for the instance */
	lb = lower(d, n, w, W);

	/* compute the TS solution */
	ub = TSpack(d, n, w, W, lb, timeLimit, &ub0, x, b, uheur);
	if (ub <= 0) {
	   printf("\n an error occurred in procedure TSpack!\n");
	   exit(0);
	} else {
		printf("\n finished packing!\n");

		for (i = 0; i < n; i++){
			printf("Piece %2d: bin: %2d, x=%2d, y=%2d\n", i, b[i], x[0][i], x[1][i]);
		}
		for (i = 0; i < ub; i++) {
			printf("\nBIN: %02x\n", i);
			for (j = 0; j < n; j++) {
				if (b[j] == i) {
  					printf("  piece %2d:, x=%2d, y=%2d, l=%2d, h=%2d\n", j, x[0][j], x[1][j], w[0][j], w[1][j]);
				}
			}
			printf("\n");
		}
	}

	/* check the correctnes of the solution */
	correct = checkfs(d, n, w, W, x, b);
	if (!correct) {
	   printf("\n the final solution is not feasible!\n");
	   exit(0);
	}
	
	/* print the solution value */
	printf(" LB = %4d iUB = %4d fUB = %4d\n", lb, ub0, ub);

	/* memory de-allocation */
	for (i = 0; i < d; i++) {
	    free(w[i]); free(x[i]);
	}
	free(w); free(x);
	free(W); free(b);

	return(0);
}
