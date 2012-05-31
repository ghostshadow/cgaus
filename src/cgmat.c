#include "cgmat.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define REPOOM(x) {fprintf(stderr,"\n\nOut of memory: %s\n\n",x); abort();}
#define REPOOMN(x,y) {fprintf(stderr,"\n\nOut of memory: %s; []=%d\n\n",x,y); abort();}

struct mat {
	double** val;
	int col,row;
};

/*
 * returns col
 */
int mat_get_col(mxp mat) {return mat->col;}

/*
 * returns row
 */
int mat_get_row(mxp mat) {return mat->row;}

/*
 * multiplies every number in a line of the length col with fact
 */
static void mat_mull(double* line,int col,double fact) {
	for(int ct=0;ct<col;ct++) line[ct]*=fact;
}

/*
 * adds to each number in a line of the length col the specific number in the second equaly length line
 */
static void mat_addl(double* tline,int col,double* sline) {
	for(int ct=0;ct<col;ct++) tline[ct]+=sline[ct];
}

/*
 * swaps the position of the current line with a later line, 
 * that is non zero at the spos position and returns 1, 
 * if no matching line is found return 0
 */
static int mat_swapl(mxp mat,int spos) {
	if((mat->val)[spos-1][spos-1]!=0) ; /*err?!*/
	double* tmpp;
	if(!(tmpp=(mat->val)[spos-1]))
		REPOOM("tmpp@mat_swapl()")
	for(int rt=spos;rt<(mat->row);rt++) {
		if((mat->val)[rt][spos-1]==0)
			continue;
		(mat->val)[spos-1]=(mat->val)[rt];
		(mat->val)[rt]=tmpp;
		return 1;
	}
	return 0;
}

/*
 * solving the matrix to format where there is only a diagonal line of ones 
 * but the rest zeros expect last column
 */
static mxp mat_solve_norm(mxp mat) {
	for(int rt=0;rt<((mat->col)-1<=(mat->row)?(mat->col)-1:(mat->row));rt++) {
		if((mat->val)[rt][rt]==0) 
			if(!mat_swapl(mat,rt+1)) 
				continue;
		mat_mull((mat->val)[rt],mat->col,-1/(mat->val)[rt][rt]);
		double* tmp;
		if(!(tmp=malloc((mat->col)*sizeof(double))))
			REPOOM("tmp(1)@mat_solve_norm()");
		for(int rt2=rt+1;rt2<(mat->row);rt2++) {
			memcpy(tmp,(mat->val)[rt],(mat->col)*sizeof(double));
			mat_mull(tmp,(mat->col),(mat->val)[rt2][rt]);
			mat_addl((mat->val)[rt2],(mat->col),tmp);
		}
		free(tmp);
	}
	for(int tr=((mat->col)-2<=(mat->row)-1?(mat->col)-2:(mat->row)-1);tr>=0;tr--) {
		if((mat->val)[tr][tr]==0) 
			continue;
		double* tmp;
		if(!(tmp=malloc((mat->col)*sizeof(double)))) 
			REPOOM("tmp(2)@mat_solve_norm()");
		for(int tr2=tr-1;tr2>=0;tr2--) {
			memcpy(tmp,(mat->val)[tr],(mat->col)*sizeof(double));
			mat_mull(tmp,(mat->col),(mat->val)[tr2][tr]);
			mat_addl((mat->val)[tr2],(mat->col),tmp);
		}
		free(tmp);
	}
	return mat;
}

/*
 * returns 1 if all lines with only zeros in the front ends with zero
 * else return 0
 */
static int mat_solvable(mxp mat) {
	int lok=1;
	for(int rt=(mat->row)-1;rt>=1;rt--) {
		if((mat->val)[rt-1][(mat->col)-1]==0) 
			continue;
		lok=1;
		for(int ct=0;ct<(mat->col)-2;ct++) {
			if((mat->val)[rt-1][ct]!=0) 
				continue;
			lok=0;
		}
		if(!lok) 
			return 0;
	}
	return 1;
}

/*
 * sets the matrix to minimum size of one row and two columns
 */
mxp mat_init() {
	mxp mat;
	if(!(mat=malloc(sizeof(struct mat))))
		REPOOM("mat@mat_init()");
	(mat->col)=2;
	(mat->row)=1;
	if(!((mat->val)=malloc(sizeof(double*))))
		REPOOM("mat->val@mat_init()");
	if(!((mat->val)[0]=malloc(2*sizeof(double))))
		REPOOM("mat->val[0]@mat_init()");
	(mat->val)[0][0]=0;
	(mat->val)[0][1]=0;
	return mat;
}

/*
 * destructs a matrix
 */
void mat_destr(mxp mat) {
	if(!(mat==NULL)) {
		for(int rt=0;rt<(mat->row);rt++) 
			free((mat->val)[rt]);
		free(mat->val);
		(mat->col)=0;
		(mat->row)=0;
		free(mat);
	}
}

/*
 * resizes the matrix to given size
 * -1 means no change here
 * warning: lastest numbers will be lost at srink
 */
mxp mat_resize(mxp mat,int nrow,int ncol) {
	if(!(ncol==-1 || ncol==(mat->col))) {
		double* tmp;
		if(!(tmp=malloc(ncol*sizeof(double))))
			REPOOM("tmp(1)@mat_resize()");
		for(int rt=0;rt<(mat->row);rt++) {
			memcpy(tmp,(mat->val)[rt],((mat->col)>ncol?ncol:(mat->col))*sizeof(double));
			free((mat->val)[rt]);
			if(!((mat->val)[rt]=malloc(ncol*sizeof(double)))) 
				REPOOMN("mat->val[rt]@mat_resize()",rt);
			memcpy((mat->val)[rt],tmp,ncol*sizeof(double));
		}
		(mat->col)=ncol;
		free(tmp);
	}
	if(!(nrow==-1 || nrow==(mat->row))) {
		double** tmp;
		if(!(tmp=malloc(nrow*sizeof(double*))))
			REPOOM("tmp(2)@mat_resize()");
		if(nrow<(mat->row))
			for(int rt=nrow;rt<(mat->row);rt++)
				free((mat->val)[rt]);
		memcpy(tmp,(mat->val),((mat->row)>nrow?nrow:(mat->row))*sizeof(double*));
		free(mat->val);
		if(!((mat->val)=malloc(nrow*sizeof(double*))))
			REPOOM("mat->val@mat_resize()");
		memcpy((mat->val),tmp,nrow*sizeof(double*));
		if(nrow>(mat->row))
			for(int rt=(mat->row);rt<nrow;rt++) 
				if(!((mat->val)[rt]=malloc((mat->col)*sizeof(double))))
					REPOOMN("mat->val[rt]@mat_resize()",rt);
		(mat->row)=nrow;
		free(tmp);
	}
	return mat;
}

/*
 * returns value at given position or 1 if fail
 */
double mat_get_val(mxp mat,int r,int c) {
	if(r<(mat->row) && c<(mat->col) && r>=0 && c>=0)
		return (mat->val)[r][c];
	return 1;
}

/*
 * sets the value of the given position to given number
 */
void mat_set_val(mxp mat,int r,int c,double val) {
	if(r<(mat->row) && c<(mat->col) && r>=0 && c>=0)
		(mat->val)[r][c]=val;
}

/*
 * starts solving with gauss methode
 * returns 1 if the system was solvable
 * or 0 when not
 */
int mat_solve(mxp mat) {
	mat_solve_norm(mat);
	return mat_solvable(mat);
}

