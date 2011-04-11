#include "TSpack.h"

/** Main entry point. */
int main(int argc, char **argv)
{                                           /* driver program for the TSPACK code */
	int    i, j, value, correct;
	int    d, n, lb, ub0, ub;
	float  timeLimit;
	char   fname[80];
	char   resultsFile[80];
	int    uh;
	FILE   *file;

	int    **w, **x;
	int    *W, *b;
	char   **id;
	char   *textline;

	if (argc < 2) {
	   printf("\n Usage: %s fileName [rotation] [outfile] [timelimit]\n",argv[0]);
	   printf("\n\n");
	   printf("File structure:\n");
	   printf("===============\n");
	   printf("d n     \t\t: d = 2[3] for x,y[,z], n = no of items\n");
	   printf("X Y [Z] \t\t: X,Y[,Z] maximum dimensions of sheet\n");
	   printf("x y [z] \t\t: n lines of x,y[,z] dims of items\n\n");
	   exit(0);
	}

	sscanf(argv[1], "%s", fname);
    timeLimit = 300;
	if (argc < 3) {
        printf("\n Allowing part rotation\n");
    } else {
	    sscanf(argv[2], "%d", &uh);
	    if (argc > 3) {
            sscanf(argv[3], "%s", resultsFile);
            if (argc > 4) {
               sscanf(argv[4], "%d", &timeLimit);
            }       
        }
    }

	file = fopen(fname, "r");
	
	if (fgetc(file) == '[') {
        // it has headers
        printf("\n HEADERS not yet implemented!\n\n");
        fclose(file);
        exit(0);
    } else {
        rewind(file);
        
        // it is plain
        fscanf(file, "%d", &d);  // dimensionality
        fscanf(file, "%d", &n);  // number of pieces

    	/* memory allocation */
    	W = (int*)calloc(d, sizeof(int));
    	b = (int*)calloc(n, sizeof(int));
    	
    	w = (int**)calloc(d, sizeof(int*));
    	x = (int**)calloc(d, sizeof(int*));
    	for (i = 0; i < d; i++) {
    	    w[i] = (int*)calloc(n, sizeof(int));
    	    x[i] = (int*)calloc(n, sizeof(int));
    	}
    	
        id = (char**)calloc(n, sizeof(char*));
        textline = (char*)calloc(121, sizeof(char));
        
        for (i = 0; i < n; i++)
            id[i] = (char*)calloc(21, sizeof(char));
          	
    
    	/* get the bin dimensions */
    	for (i = 0; i < d; i++) {
    	    fscanf(file, "%d", &value);
    	    W[i] = value;
    	}
    	/* get the item dimensions */
    	for (i = 0; i < n; i++) {
    	    for (j = 0; j < d; j++) {
    	        fscanf(file, "%d", &value);
    	        w[j][i] = value;
    	    }
    	    
    	    /* get item descriptor, if any */
    	    textline = fgets(textline, 120, file);
    	    for (j = 0; j < 20; j++) {
                if (textline[j] == '\n') {
                      j = 20;
                } else {
                   id[i][j] = textline[j];
                }
            }
        }
        fclose(file);
    }
	
    printf("\n Solving %dD problem with %d pieces...\n",d, n);

	/* compute an initial lower bound for the instance */
	lb = lower(d, n, w, W);
    printf(" OK, lower bound is %d bins of %d x %d\n", lb, W[0], W[1]);

	/* compute the TS solution */
	ub = TSpack(d, n, w, W, lb, timeLimit, &ub0, x, b, uh);
	if (ub <= 0) {
	   printf("\n an error occurred in procedure TSpack!\n");
	   exit(0);
	} else
	  printf("\n Done, now checking results... \n");

	/* check the correctnes of the solution */
	correct = checkfs(d, n, w, W, x, b);
	if (!correct) {
	   printf("\n the final solution is not feasible!\n");
	   exit(0);
	}

	/* print the solution value */
	printf(" LB = %4d iUB = %4d fUB = %4d\n", lb, ub0, ub);
    printf("\n\n---\n%d sheets required\n", ub);
    
    if (argc > 3) {
        file = fopen(resultsFile, "w");
        fprintf(file, "%d sheets of %d x %d required\n", ub, W[0], W[1]);
        for (i = 0; i < ub; i++) {
            fprintf(file, "\n\nLAYOUT: %03d\n===========\n", i + 1);
            for (j = 0; j < n; j++) {
                if (b[j] == i) fprintf(file, "  part: %6d, x: %6d, y: %6d, px: %6d, py: %6d, id: %s\n", 
                                          j + 1, x[0][j], x[1][j], w[0][j], w[1][j], id[j]);
            }
        }
        fclose(file);
    } else {
        file = fopen(fname, "a");
    	if (file != NULL) {
            fprintf( file, "\n\n---\n%d sheets required\n", ub);
            fclose(file);
        }
    }

	/* memory de-allocation */
	for (i = 0; i < d; i++) {
	    free(w[i]); free(x[i]);
	}
	free(w); free(x);
	free(W); free(b);

	return(0);
}
