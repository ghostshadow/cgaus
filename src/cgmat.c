#include "cgmat.h"

struct mat {
	double** val;
	int col,row;
};

/*
 * returns col
 */
int mat_get_col(void* mat) {return ((struct mat*)mat)->col;}

/*
 * returns row
 */
int mat_get_row(void* mat) {return ((struct mat*)mat)->row;}

/*
 * multiplies every number in a line of the length col with fact
 */
static void mat_mull(double* line,int col,double fact) {
	for(int ct=0;ct<col;ct++) *(line+ct)*=fact;
}

/*
 * adds to each number in a line of the length col the specific number in the second equaly length line
 */
static void mat_addl(double* tline,int col,double* sline) {
	for(int ct=0;ct<col;ct++) *(tline+ct)+=*(sline+ct);
}

/*
 * swaps the position of the current line with a later line, 
 * that is non zero at the spos position and returns 1, 
 * if no matching line is found return 0
 */
static int mat_swapl(void* mat,int spos) {
	if(*(*((((struct mat*)mat)->val)+spos-1)+spos-1)!=0) ; /*err?!*/
	double* tmpp=*((((struct mat*)mat)->val)+spos-1);
	for(int rt=spos;rt<(((struct mat*)mat)->row);rt++) {
		if(*(*((((struct mat*)mat)->val)+rt)+spos-1)==0)
			continue;
		*((((struct mat*)mat)->val)+spos-1)=*((((struct mat*)mat)->val)+rt);
		*((((struct mat*)mat)->val)+rt)=tmpp;
		return 1;
	}
	return 0;
}

/*
 * solving the matrix to format where there is only a diagonal line of ones 
 * but the rest zeros expect last column
 */
static void* mat_solve_norm(void* mat) {
	for(int rt=0;rt<(((struct mat*)mat)->col)-1;rt++) {
		if(*(*((((struct mat*)mat)->val)+rt)+rt)==0) 
			if(!mat_swapl(mat,rt+1)) 
				continue;
		mat_mull(*((((struct mat*)mat)->val)+rt),(((struct mat*)mat)->col),-1/(*(*((((struct mat*)mat)->val)+rt)+rt)));
		double* tmp=malloc((((struct mat*)mat)->col)*sizeof(double));
		for(int rt2=rt+1;rt2<(((struct mat*)mat)->row);rt2++) {
			memcpy(tmp,*((((struct mat*)mat)->val)+rt),(((struct mat*)mat)->col));
			mat_mull(tmp,(((struct mat*)mat)->col),*(*((((struct mat*)mat)->val)+rt2)+rt));
			mat_addl(*((((struct mat*)mat)->val)+rt2),(((struct mat*)mat)->col),tmp);
		}
		free(tmp);
	}
	for(int tr=(((struct mat*)mat)->col)-1;tr>=0;tr--) {
		if(*(*((((struct mat*)mat)->val)+tr)+tr)==0) 
			continue;
		double* tmp=malloc((((struct mat*)mat)->col)*sizeof(double));
		for(int tr2=tr-1;tr2>=0;tr2--) {
			memcpy(tmp,*((((struct mat*)mat)->val)+tr),(((struct mat*)mat)->col));
			mat_mull(tmp,(((struct mat*)mat)->col),*(*((((struct mat*)mat)->val)+tr2)+tr));
			mat_addl(*((((struct mat*)mat)->val)+tr2),(((struct mat*)mat)->col),tmp);
		}
		free(tmp);
	}
	return mat;
}

/*
 * returns 1 if all lines with only zeros in the front ends with zero
 * else return 0
 */
static int mat_solvable(void* mat) {
	for(int rt=(((struct mat*)mat)->row);rt>=1;rt--) {
		if(*(*((((struct mat*)mat)->val)+rt-1)+(((struct mat*)mat)->col)-1)==0) 
			continue;
		for(int ct=0;ct<(((struct mat*)mat)->col)-1;ct++) 
			if(*(*((((struct mat*)mat)->val)+rt-1)+ct)!=0) 
				return 0;
	}
	return 1;
}

/*
 * sets the matrix to minimum size of one row and two columns
 */
void* mat_init(void* mat) {
	(((struct mat*)mat)->col)=2;
	(((struct mat*)mat)->row)=1;
	(((struct mat*)mat)->val)=malloc(sizeof(double*));
	*(((struct mat*)mat)->val)=malloc(2*sizeof(double));
	return mat;
}

/*
 * destructs a matrix
 */
void mat_destr(void* mat) {
	if(!(mat==NULL)) {
		for(int rt=0;rt<(((struct mat*)mat)->row);rt++) 
			free(*((((struct mat*)mat)->val)+rt));
		free(((struct mat*)mat)->val);
		(((struct mat*)mat)->col)=0;
		(((struct mat*)mat)->row)=0;
		mat=NULL;
	}
}

/*
 * resizes the matrix to given size
 * -1 means no change here
 * warning: lastest numbers will be lost at srink
 */
void* mat_resize(void* mat,int ncol,int nrow) {
	if(!(ncol==-1 || ncol==(((struct mat*)mat)->col))) {
		double* tmp=malloc(ncol*sizeof(double));
		for(int rt=0;rt<(((struct mat*)mat)->row);rt++) {
			memcpy(tmp,*((((struct mat*)mat)->val)+rt),((((struct mat*)mat)->col)>ncol?ncol:(((struct mat*)mat)->col)));
			*((((struct mat*)mat)->val)+rt)=realloc(*((((struct mat*)mat)->val)+rt),ncol*sizeof(double));
			memcpy(*((((struct mat*)mat)->val)+rt),tmp,ncol);
		}
		(((struct mat*)mat)->col)=ncol;
		free(tmp);
	}
	if(!(nrow==-1 || nrow==(((struct mat*)mat)->row))) {
		double** tmp=malloc(nrow*sizeof(double*));
		memcpy(tmp,(((struct mat*)mat)->val),((((struct mat*)mat)->row)>nrow?nrow:(((struct mat*)mat)->row)));
		(((struct mat*)mat)->val)=realloc((((struct mat*)mat)->val),nrow*sizeof(double*));
		memcpy((((struct mat*)mat)->val),tmp,nrow);
		if(nrow>(((struct mat*)mat)->row))
			for(int rt=(((struct mat*)mat)->row);rt<nrow;rt++) 
				*((((struct mat*)mat)->val)+rt)=malloc((((struct mat*)mat)->col)*sizeof(double));
		(((struct mat*)mat)->row)=nrow;
		free(tmp);
	}
	return mat;
}

/*
 * returns value at given position
 */
double mat_get_val(void* mat,int r,int c) {return *(*((((struct mat*)mat)->val)+r-1)+c-1);}

/*
 * sets the value of the given position to given number
 */
void mat_set_val(void* mat,int r,int c,double val) {*(*((((struct mat*)mat)->val)+r-1)+c-1)=val;}

/*
 * starts solving with gauss methode
 * returns 1 if the system was solvable
 * or 0 when not
 */
int mat_solve(void* mat) {
	mat_solve_norm(mat);
	return mat_solvable(mat);
}

