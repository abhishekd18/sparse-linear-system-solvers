#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mmio.h"
#include "functions.h"

/*
*	GMRES implementation
*	
*	MAT		:	Data structure of type MTX
*	x0 		:	Initial guess vector
*	xm 		:	Computed Solution vector
*	b 		:	Right hand side vector
*	m 		:	No. of Krylov vectors/No. of iterations
*	tol		:	Stopping tolerence
*	res 		:	Type of residual to be monitored for convergence- "Absolute" or "Relative"
*	preconditioner	: 	Type of preconditioner to be applied
*	mode 		:	Operation mode for GMRES- "full" or "restarted"
*/

double GMRES(MTX *MAT, double* x0, double* xm, double* b, int *m, double tol, char* res, char* preconditioner, char* mode){
	
	if(strcmp(mode,"full")==0)	*m = 10000;
	
	// Allocate memory for residual vector r0
	double* r0 = (double*) malloc(MAT->ncols*sizeof(double));

	// Allocate memory for g
	double* g = (double*) malloc((*m+1)*sizeof(double));

	// Allocate memory for V
	double** V = (double**) malloc((*m+1)*sizeof(double*));
	for(int j=0;j<*m+1;j++)
		V[j] = (double*) malloc(MAT->ncols*sizeof(double));

	// Allocate memory for Hm
	double** H = AllocateDynamicArray(*m+1,*m);

	// Allocate memory for Rm
	double **R = AllocateDynamicArray(*m+1,*m);

	double Rk, Rk_1;
	double *c = (double*) malloc(*m*sizeof(double));
	double *s = (double*) malloc(*m*sizeof(double));

	double *error = (double*) malloc(MAT->ncols*sizeof(double));
	double rho = 1.0, Err, rresidual, residual, residual0;
	int flag = 0;

	// Calculate Initial residual
	//r0 = A*x0
	Mat_Vec_Mult(MAT,x0,r0);
	
	//r0 = b - A*x0
	for(int i=0;i<MAT->ncols;i++)	
		r0[i] = b[i] - r0[i];

	// res = ||r0||
	residual0 = vecnorm(r0,MAT->ncols);

	// residual for internal for loop
	residual = residual0;

	FILE *fp;
	fp = fopen("./Output/Error_GMRES.out","a");

	while(rho>tol){

		// Preconditioning
		if(strcmp(preconditioner,"NULL")!=0){
			LeftPreconditioning(MAT,r0,preconditioner);
			residual = vecnorm(r0,MAT->ncols);
		}		
	
		g[0] = vecnorm(r0,MAT->ncols);
		
		for(int i=1;i<*m+1;i++)
			g[i] = 0;

		// Calculate V[1] as r0/||r0||
		for(int i=0;i<MAT->ncols;i++)
			V[0][i] = r0[i]/residual;

		for(int j=0;j<*m;j++){

			get_Krylov(MAT,V,H,j,preconditioner);

			R[0][j] = H[0][j];
			for(int k=1;k<=j;k++){
				Rk_1 = c[k-1]*R[k-1][j] + s[k-1]*H[k][j];
				Rk = -s[k-1]*R[k-1][j] + c[k-1]*H[k][j];
			
				R[k-1][j] = Rk_1;
				R[k][j] = Rk;
			}
			c[j] = R[j][j]/pow((R[j][j]*R[j][j] + H[j+1][j]*H[j+1][j]),0.5);
			s[j] = H[j+1][j]/pow((R[j][j]*R[j][j] + H[j+1][j]*H[j+1][j]),0.5);
			R[j][j] = c[j]*R[j][j] + s[j]*H[j+1][j];	
	
			g[j+1] = -s[j]*g[j];
			g[j] = c[j]*g[j];

			if(strcmp(res,"Relative")==0){
				rresidual = fabs(g[j+1])/residual;
				if(strcmp(mode,"full")==0)
					fprintf(stdout,"m = %d\tRelative Residual = %lg\n",j+1,rresidual);
				else
					fprintf(stdout,"m = %d\tRelative Residual = %lg\n",j+1,rresidual*residual/residual0);
			}else if(strcmp(res,"Absolute")==0){
				rresidual = fabs(g[j+1]);
				fprintf(stdout,"m = %d\tAbsolute Residual = %lg\n",j+1,rresidual);
			}else{
				fprintf(stderr,"\nPlease check the type of residual to be specified. \
					Keywords to be used are:\n \"Absolute\"\tOR\t\"Relative\"");
				exit(1);
			}

			if(strcmp(mode,"full")==0){
				if(rresidual<tol){	
					*m=j;	break;
				}
			}	

			// Calculate xm = x0 + V*y
			get_solution(MAT, R, V, g, j, x0, xm);

			// error = x - xm
			for(int i=0;i<MAT->ncols;i++){
				error[i] = 1.0 - xm[i];
			}

			// Err = ||x-xm||
			Err = vecnorm(error,MAT->ncols);

			// Write :	Iteration No.	Error-2-Norm	Rel Residual	Abs Residual
			if(strcmp(res,"Relative")==0){
				if(strcmp(mode,"full")==0)
				fprintf(fp,"%d\t%1.16E\t%1.16E\t%1.16E\n", (*m*flag+j)+1, Err, rresidual, rresidual*residual);
				else
				fprintf(fp,"%d\t%1.16E\t%1.16E\t%1.16E\n", (*m*flag+j)+1, Err, rresidual*residual/residual0, rresidual*residual);
			}else if(strcmp(res,"Absolute")==0){
				fprintf(fp,"%d\t%1.16E\t%1.16E\t%1.16E\n", (*m*flag+j)+1, Err, rresidual/residual0, rresidual);
			}

		}// for end


		if(strcmp(res,"Relative")==0){
			if(strcmp(mode,"full")==0)
				rho = rresidual;
			else			
				rho = rresidual*residual/residual0;		
		}else if(strcmp(res,"Absolute")==0){
			rho = rresidual;
		}

		flag++;

		// Calculate residual for next iteration
		for(int i=0;i<MAT->ncols;i++)
			x0[i] = xm[i];

		//r = A*x0
		Mat_Vec_Mult(MAT,x0,r0);
	
		//r = b - A*x0
		for(int i=0;i<MAT->ncols;i++)	
			r0[i] = b[i] - r0[i];

		// res = ||r0||
		residual = vecnorm(r0,MAT->ncols);

	}// While end

	fclose(fp);

	double eps = 1.0e-3;	// Tolerance check
	check_orthonormality(V, *m, MAT->nrows, eps);
	
	FreeDynamicArray(R,*m+1);
	FreeDynamicArray(H,*m+1);
	
	for(int j=0;j<*m+1;j++)
		free(V[j]);
	free(V);

	// Assign cumulative value to m in case of restarted mode. 
	int iter = (*m+1)*flag;
	*m = iter;

	free(error);
	free(r0);
	free(c);
	free(s);
	free(g);

return rho;
}

/*
*	Get solution vector
*/
void get_solution(MTX *MAT, double** R, double** V, double* g, int j, double* x0, double* xm){

	// Allocate memory for y
	double* y = (double*) malloc((j+1)*sizeof(double));
	double* sum = (double*) malloc(MAT->ncols*sizeof(double));

	// y = Inv(Rm)*g
	Back_Substitute(R, j+1, j+1, g, y); 


	for(int i=0;i<MAT->ncols;i++)
		sum[i] = 0.0;

	// V*y
	for(int i=0;i<MAT->ncols;i++)
		for(int k=0;k<j+1;k++)
			sum[i] = sum[i]  + V[k][i]*y[k];

	// xm = x0 + V*y
	for(int i=0;i<MAT->ncols;i++)
		xm[i] = x0[i] + sum[i];

	free(y);
	free(sum);
return;
}

/*
*	Performs Gram-Schmidt algorithm and gives out orthonormal Krylov vectors
*/
void get_Krylov(MTX *MAT, double** V, double** H, int j, char* preconditioner){

	double* w = (double*) malloc(MAT->ncols*sizeof(double));

	//w=A*vj
	Mat_Vec_Mult(MAT,V[j],w);

	// Preconditioning w=Inv(M)*w
	if(strcmp(preconditioner,"NULL")!=0)	LeftPreconditioning(MAT,w,preconditioner);

	for(int i=0;i<=j;i++){
		H[i][j]=scalarProd(V[i],w,MAT->ncols);	
		for(int k=0;k<MAT->ncols;k++)
			w[k] = w[k] - H[i][j]*V[i][k];	
	}

	H[j+1][j] = vecnorm(w,MAT->ncols);

	for(int k=0;k<MAT->ncols;k++)
		V[j+1][k] = w[k]/H[j+1][j];

	free(w);
return;
}

/*
*	Finds the norm of a vector of size n
*/
double vecnorm(double* vec, int n){
	double norm = 0.0;
	for(int i=0;i<n;i++)
		norm = norm + vec[i]*vec[i];
	norm = pow(norm,0.5);
return norm;
}

/*
*	Finds the scalar product of vec1 and vec2 of size n
*/
double scalarProd(double* vec1, double* vec2, int n){
	double prod = 0.0;
	for(int i=0;i<n;i++)
		prod = prod + vec1[i]*vec2[i];
return prod;
}

/*
*	Back substitution function for system with upper triangular matrix
*/
void Back_Substitute(double** R, int nRows, int nCols, double* invec, double* outvec){
	
	double* temp = (double*) malloc(nRows*sizeof(double)); 
	for(int i=0;i<nRows;i++)
		temp[i] = invec[i];

	outvec[nRows-1] = temp[nRows-1]/R[nRows-1][nCols-1];
	for(int i=nRows-2;i>=0;i--){
		for(int j=nCols-1;j>i;j--)
			temp[i] = temp[i] - R[i][j]*outvec[j];
		
		outvec[i] = temp[i]/R[i][i];
	}
	
	free(temp);
return;
}
