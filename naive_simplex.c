/* 
   The naive simplex algorithm as presented in Ch. 2
   of our textbook

   Usage: 
   To compile: gcc -O3 naive_simplex.c -o simplex
   Then: simplex  input_file_name 
   or:  simplex -v input_file_name 
   for more detailed output.

   Input file format: a sequence of numbers separated 
   by space,  tab or end-of-line, starting with the number 
   of variables n, the number of constrains m, then the 
   augmented matrix of the equation constraints [A b], 
   then the coefficients of the objective function that 
   is being maximized.
   
   Note that this code is just an illustration and has not 
   been optimized for speed, and I haven't learned a lot 
   of CS so the code quality is sloppy :)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int verbose=0;

typedef struct
{
  int m, n;
  double* d;
} aug_mat;

int print_aug_mat(aug_mat M){
  int i,j;
  for(i=0;i<M.m;i++){
    for(j=0;j<M.n;j++)
      printf("%g ",M.d[i*M.n+j]);
    printf("\n");
  }
  return 0;
}

/* Permute the rows i and j of matrix M */
int perm(int i, int j, aug_mat M){
  double buf;
  double* c=M.d;
  if(i==j)
    return 0;
  int k;
  for(k=0;k<M.n;k++){
    buf=c[i*M.n+k];
    c[i*M.n+k]=c[j*M.n+k];
    c[j*M.n+k]=buf;
  }
  return 0;
}

/* Turn the i-th column of the matrix into standard 
 vector e_k via row operations*/
int elim(int i, int k, aug_mat M){
  double* c=M.d;
  int j;
  for(j=0;j<M.m;j++){
    if(j!=k){
      double t=c[j*M.n+i];
      /*printf("%d,%d\n",j,k);*/
      int s;
      for(s=0;s<M.n;s++)
	c[j*M.n+s]-=c[k*M.n+s]*t/c[k*M.n+i];
    }
  }
  double t=c[k*M.n+i];
  for(j=0;j<M.n;j++)
	c[k*M.n+j]/=t;
    
  return 0;
}


int is_basic(int i, int m, int* basics){
  int j;
  for(j=0;j<m;j++)
    if(i==basics[j])return 1;
  return 0;
}

/* Simplex algorithm starting from basic variables */
int simplex_from_basic(int* basics,  double* o_coeffs,
		       aug_mat M){
  int n=M.n-1;
  int m=M.m;
  double* c=M.d;
  double* cur_o_coeffs=(double*)malloc(n*sizeof(double));
  double* old_o_coeffs=(double*)malloc(n*sizeof(double));
  /*Initialize */
  int i,j;
  for(i=0;i<n;i++){
    if(is_basic(i,m,basics))cur_o_coeffs[i]=0;
    else{
      cur_o_coeffs[i]=o_coeffs[i];
      for(j=0;j<m;j++)
	cur_o_coeffs[i]-=o_coeffs[basics[j]]*c[j*M.n+i];
    }
  }

  /*Change basic variables */
  while(1){
    if(verbose){
      for(i=0;i<m;i++)printf("%d,",basics[i]);
      printf("\n");
      print_aug_mat(M);
      for(i=0;i<n;i++)printf("%g ",cur_o_coeffs[i]);
      printf("\n\n");
    }
    for(i=0;i<n;i++)
      if(cur_o_coeffs[i]>0)break;
    if(i==n)
      return 0;
    for(j=0;j<m;j++)
      if(c[j*M.n+n]*c[j*M.n+i]>0)break;
    if(j==m)
      return 1;
    double opt_ratio=c[j*M.n+n]/c[j*M.n+i];
    int outv=j;
    for(;j<m;j++)
      if(c[j*M.n+n]*c[j*M.n+i]>0&&
	 opt_ratio>c[j*M.n+n]/c[j*M.n+i]){
	outv=j;
	opt_ratio=c[j*M.n+n]/c[j*M.n+i];
      }
    elim(i,outv,M);
    int k;
    for(k=0;k<n;k++)
      old_o_coeffs[k]=cur_o_coeffs[k];
    for(k=0;k<n;k++)
      cur_o_coeffs[k]-=old_o_coeffs[i]*c[outv*M.n+k];
    cur_o_coeffs[i]=0;
    basics[outv]=i;
  }
  free(cur_o_coeffs);
  free(old_o_coeffs);
  return 0;
}

int main(int argc, char* argv[]){

  /* Read data */
  char* input_filename=argv[argc-1];
  if(argv[1][0]=='-'&&argv[1][1]=='v'){
    verbose=1;
  }
  FILE* data=fopen(input_filename, "r");
  int m, n, i, j;
  if(!fscanf(data, "%d", &n))return 1;
  if(!fscanf(data, "%d", &m))return 1;
  double* c_coeffs=(double*)malloc(m*(n+m+1)*sizeof(double));
  double* o_coeffs=(double*)malloc(n*sizeof(double));
  for(i=0;i<m;i++){
    for(j=0;j<n;j++)
      if(!fscanf(data,"%lf",c_coeffs+i*(n+m+1)+j))return 1;
    if(!fscanf(data,"%lf",c_coeffs+i*(n+m+1)+m+n))return 1;
    double b=c_coeffs[i*(n+m+1)+m+n];
    for(j=n;j<m+n;j++)
      if(j-n==i){
	if(b>=0)c_coeffs[i*(n+m+1)+j]=1;
	else c_coeffs[i*(n+m+1)+j]=-1;
      }
      else c_coeffs[i*(n+m+1)+j]=0;
  }
    
  for(i=0;i<n;i++)
    if(!fscanf(data,"%lf",o_coeffs+i))return 1;
  fclose(data);


  int* basics=(int*)malloc(m*sizeof(int));
  /*phase 1 */
  if(verbose)printf("Phase 1\n");
  double* o_coeffs1=(double*)malloc((n+m)*sizeof(double));
  for(i=0;i<n;i++)
    o_coeffs1[i]=0;
  for(i=n;i<m+n;i++)
    o_coeffs1[i]=-1;
  for(i=0;i<m;i++)
    basics[i]=n+i;
  aug_mat M;
  M.m=m;
  M.n=m+n+1;
  M.d=c_coeffs;
  simplex_from_basic(basics, o_coeffs1, M);
  for(i=0;i<m;i++)
    if(basics[i]>=n){
      printf("No feasible solution.\n");
      break;
    }

  /*phase 2*/
  if(i==m){
    if(verbose)printf("Phase 2\n");
    double* c_coeffs2=(double*)malloc(m*(n+1)*sizeof(double));
    for(i=0;i<m;i++){
      for(j=0;j<n;j++)
	c_coeffs2[i*(n+1)+j]=c_coeffs[i*(m+n+1)+j];
      c_coeffs2[i*(n+1)+n]=c_coeffs[i*(m+n+1)+m+n];
    }
    free(c_coeffs);
    M.d=c_coeffs2;
    M.n=n+1;
    int r=simplex_from_basic(basics, o_coeffs, M);
    if(r==1)printf("No optimal solution.\n");
    else {
      double* optv=o_coeffs1;
      double opt=0;
      for(i=0;i<n;i++)optv[i]=0;
      for(j=0;j<m;j++){
	optv[basics[j]]=M.d[j*(n+1)+n];
        opt+=M.d[j*(n+1)+n]*o_coeffs[basics[j]];
      }
      printf("Optimal value: %g \n Optimal solution:\n",opt);
      for(i=0;i<n;i++)
	printf("%g\n",optv[i]);
    }
  }

  return 0;
}
