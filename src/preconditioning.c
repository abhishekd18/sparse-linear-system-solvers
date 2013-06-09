#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mmio.h"
#include "functions.h"

/*
*	Left Preconditioner
*/
void LeftPreconditioning(MTX *MAT, double* vec, char* preconditioner){

	if(strcmp(preconditioner,"JACOBI") == 0){

		for(int j=0;j<MAT->nrows;j++)
			for(int i=MAT->row_ptr[j];i<MAT->row_ptr[j+1];i++)
				if(MAT->JA[i]==j)
					vec[j] = vec[j]/MAT->val[i];
	
	}else if(strcmp(preconditioner,"GAUSS_SEIDEL") == 0){

		int i;
		for(int j=0;j<MAT->nrows;j++){
			for(i=MAT->row_ptr[j];i<MAT->row_ptr[j+1];i++){
				if(MAT->JA[i]==j)	break;
				vec[j] = vec[j] - MAT->val[i]*vec[MAT->JA[i]];
			}
			vec[j] = vec[j]/MAT->val[i];
		}

	}else
		fprintf(stderr,"\nPlease check the choice of preconditioner again!\n\
			Available options are:\t1. JACOBI\t2. GAUSS_SEIDEL\t3. NULL\n");

return;
}
