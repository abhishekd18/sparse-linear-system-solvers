#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mmio.h"
#include "functions.h"

/*
*	Conjugate Gradient
*/
double CG(MTX *MAT, double* x0, double* xm, double* b, int *iter, double tol){
	*iter = 0;

	/* Allocate memory for residual vector rm */
	double* rm = (double*) malloc(MAT->ncols*sizeof(double));

	//r0 = A*x0
	Mat_Vec_Mult(MAT,x0,rm);
	
	//r0 = b - A*x0
	for(int i=0;i<MAT->ncols;i++)	
		rm[i] = b[i] - rm[i];

	// xm = x0 for m=1
	for(int i=0;i<MAT->ncols;i++)	
		xm[i] = x0[i];

	/* Allocate memory for residual vector pm */
	double* pm = (double*) malloc(MAT->ncols*sizeof(double));

	//p0 = r0
	for(int i=0;i<MAT->ncols;i++)	
		pm[i] = rm[i];

	double* Apm = (double*) malloc(MAT->ncols*sizeof(double));
	
	double Apm_pm, rm_rm_old, rm_rm_new, residual=1.0, residual0, alpha_m, beta_m, Err_Anorm, Err;

	//(r0,r0)
	rm_rm_old = scalarProd(rm,rm,MAT->ncols);
	residual0 = pow(rm_rm_old,0.5);
	
	double* error = (double*) malloc(MAT->ncols*sizeof(double));
	double* Aerror = (double*) malloc(MAT->ncols*sizeof(double));

	// Print Error Vs. Iteration in a file
	FILE* fp = fopen("./Output/Error_CG_Anorm.out","a");
	FILE* fp1 = fopen("./Output/Error_CG.out","a");

	while((residual/residual0)>tol){

		// Count iterations
		*iter = *iter + 1;

		//A*pm
		Mat_Vec_Mult(MAT,pm,Apm);

		//(A*pm,pm)
		Apm_pm = scalarProd(Apm,pm,MAT->ncols);

		//alpha = (rm,rm)/(A*pm,pm)
		alpha_m = rm_rm_old/Apm_pm;
		
		//x_{m+1} = x_{m} + alpha_m * pm 
		for(int j=0;j<MAT->ncols;j++)
			xm[j] = xm[j] + alpha_m*pm[j];

		//r_{m+1} = rm - alpha_m * A*pm
		for(int j=0;j<MAT->ncols;j++)
			rm[j] = rm[j] - alpha_m*Apm[j];

		//(r_m+1,r_m+1)
		rm_rm_new = scalarProd(rm,rm,MAT->ncols);
		residual = pow(rm_rm_new,0.5);
	
		fprintf(stdout,"\nIteration No. = %d\tRelative Residual = %e\n", *iter, residual/residual0);

		//beta_m = -(r_m+1,r_m+1)/(rm,rm)
		beta_m = -rm_rm_new/rm_rm_old;
		
		rm_rm_old = rm_rm_new;
		
		//p_m+1 = r_m+1 - beta_m * pm
		for(int j=0;j<MAT->ncols;j++)
			pm[j] = rm[j] - beta_m*pm[j];

		// Error (x* - xm)
		for(int j=0;j<MAT->ncols;j++)
			error[j] = 1.0 - xm[j];

		// A*error
		Mat_Vec_Mult(MAT,error,Aerror);

		// Error ||x* - xm|| in "A" norm
		Err_Anorm = scalarProd(Aerror,error,MAT->ncols);
		Err_Anorm = pow(Err_Anorm,0.5);
		Err = vecnorm(error,MAT->nrows);

		fprintf(fp,"%d\t%1.16E\n", *iter, Err_Anorm);
		fprintf(fp1,"%d\t%1.16E\t%1.16E\n", *iter, Err, residual/residual0);
	}
	
	fclose(fp);
	fclose(fp1);

	free(rm);
	free(pm);
	free(Apm);
	free(error);
	free(Aerror);
	
return (residual/residual0);
}
