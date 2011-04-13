#include "TSpack.h"

/** Main entry point. */
int main(int argc, char **argv)
{                                           /* driver program for the TSPACK code */
	int    i, j, value, correct;
	int    d, n, lb, ub0, ub;
	float  yield, piece;
	float  timeLimit;
	int    iterLimit;
	char   fname[80];
	char   resultsFile[80];
	int    uh;
	FILE   *file;

	int    **w, **x;
	int    *W, *b;
	char   **id;
	char   *textline;

	if (argc < 2) {
	   printf("\n Usage: %s fileName [rotation] [outfile] [timelimit] [iterLimit]\n",argv[0]);
	   printf("\n\n");
	   printf("File structure:\n");
	   printf("===============\n");
	   printf("d n     \t\t: d = 2[3] for x,y[,z], n = no of items\n");
	   printf("X Y [Z] \t\t: X,Y[,Z] maximum dimensions of sheet\n");
	   printf("x y [z] \t\t: n lines of x,y[,z] dims of items\n\n");
	   exit(0);
	}

	sscanf(argv[1], "%s", fname);
    iterLimit = nITmax;
	if (argc < 3) {
        uh = 1;
    } else {
	    sscanf(argv[2], "%d", &uh);
    }
    if (argc > 3) {
        sscanf(argv[3], "%s", resultsFile);
    } else {
        sscanf(argv[1], "%s", resultsFile);
        textline = strrchr(resultsFile, '.');
        if (textline) *textline = 0;
        sprintf(resultsFile, "%s.out", resultsFile);
    }
    if (argc > 4) {
       sscanf(argv[4], "%f", &timeLimit);
    } else {
       timeLimit = 120;
    }

	file = fopen(fname, "r");
	printf("\nReading data from file %s...\n", fname);
	
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
    	yield = 0;
    	for (i = 0; i < n; i++) {
            piece = 1.0;
    	    for (j = 0; j < d; j++) {
    	        fscanf(file, "%d", &value);
    	        w[j][i] = value;
    	        piece = piece * (float) w[j][i] / (float) W[j];
    	    }
    	    
    	    yield = yield + piece;
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

    if (argc > 5) {
       sscanf(argv[5], "%i", &iterLimit);
    } else {
       if (n > 50) iterLimit = 6.0 * (float) n;
    }
	
    printf("\n Solving %dD problem with %d pieces %s rotation.\n",d, n, uh ? "with" : "without");
    printf(" Up to %d attempts in %.0f seconds\n\n Press 'Q' to exit early.\n\n\n", iterLimit, timeLimit);
    

	/* compute an initial lower bound for the instance */
	lb = lower(d, n, w, W);
    printf(" Minimum is %d sheets of %d x %d (%2.2f%% yield)\n", lb, W[0], W[1], yield * 100 / (float) lb);

	/* compute the TS solution */
	ub = TSpack(d, n, w, W, lb, timeLimit, &ub0, x, b, uh, iterLimit);
	if (ub <= 0) {
	   printf("\n an error occurred in procedure TSpack!\n");
	   exit(0);
	}

	/* check the correctnes of the solution */
	correct = checkfs(d, n, w, W, x, b);
	if (!correct) {
	   printf("\n the final solution is not feasible!\n");
	   exit(0);
	}

	/* print the solution value */
    if (strlen(fname) > 3) {
        file = fopen(resultsFile, "w");
        fprintf(file, "%d sheets of %d x %d required\n", ub, W[0], W[1]);
        timeLimit = 0;
        for (i = 0; i < ub; i++) {
            fprintf(file, "\n\nLAYOUT: %03d\n===========\n", i + 1);
            yield = 0;
            value = 0;
            for (j = 0; j < n; j++) {
                if (b[j] == i) {
                   fprintf(file, "  part: %6d, x: %6d, y: %6d, px: %6d, py: %6d, id: %s\n", 
                                          j + 1, x[0][j], x[1][j], w[0][j], w[1][j], id[j]);
                   yield = yield + ((float) w[0][j]) * ((float) w[1][j]);
                   value++;
                }
            }
            yield = yield * 100 / (((float) W[0]) * ((float) W[1]));
            timeLimit += yield;
            fprintf(file, "\n %i items, YIELD: %2.2f%%\n", value, yield);
        }
        timeLimit = timeLimit / (float) ub;
        fprintf(file, "\n\n=====\nTOTAL %d pieces on %d sheets, average yield %2.2f%%\n", n, ub, timeLimit);
        fclose(file);
    }
    printf("\n\n=====\nTOTAL %d pieces on %d sheets, average yield %2.2f%%\n", n, ub, timeLimit);
    printf("Results saved to %s\n", resultsFile);
    printf("\n PRESS ANY KEY TO EXIT... \n");
    getch();

	/* memory de-allocation */
	for (i = 0; i < d; i++) {
	    free(w[i]); free(x[i]);
	}
	free(w); free(x);
	free(W); free(b);

	return(0);
}
