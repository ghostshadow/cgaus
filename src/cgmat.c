#include "cgmat.h"

struct mat {
	double** val;
	int col,row;
};

/*
 * returns col
 */
int inline mat_get_col(void* mat) {return mat->col;}

/*
 * returns row
 */
int inline mat_get_row(void* mat) {return mat->row;}

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
	if(*(*((mat->val)+spos-1)+spos-1)!=0) ; /*err?!*/
	double* tmpp=*((mat->val)+spos-1);
	for(int rt=spos;rt<(mat->row);rt++) {
		if(*(*((mat->val)+rt)+spos-1)==0)
			continue;
		*((mat->val)+spos-1)=*((mat->val)+rt);
		*((mat->val)+rt)=tmpp;
		return 1;
	}
	return 0;
}

/*
 * solving the matrix to format where there is only a diagonal line of ones 
 * but the rest zeros expect last column
 */
static void* mat_solve_norm(void* mat) {
	for(int rt=0;rt<(mat->col)-1;rt++) {
		if(*(*((mat->val)+rt)+rt)==0) 
			if(!mat_swapl(mat,rt+1)) 
				continue;
		mat_mull(*((mat->val)+rt),(mat->col),-1/(*(*((mat->val)+rt)+rt)));
		double* tmp=malloc((mat->col)*sizeof(double));
		for(int rt2=rt+1;rt2<(mat->row);rt2++) {
			memcpy(tmp,*((mat->val)+rt),(mat->col));
			mat_mull(tmp,(mat->col),*(*((mat->val)+rt2)+rt));
			mat_addl(*((mat->val)+rt2),(mat->col),tmp);
		}
		free(tmp);
	}
	for(int tr=col-1;tr>=0;tr--) {
		if(*(*((mat->val)+tr)+tr)==0) 
			continue;
		double* tmp=malloc((mat->col)*sizeof(double));
		for(int tr2=tr-1;tr2>=0;tr2--) {
			memcpy(tmp,*((mat->val)+tr),(mat->col));
			mat_mull(tmp,(mat->col),*(*((mat->val)+tr2)+tr));
			mat_addl(*((mat->val)+tr2),(mat->col),tmp);
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
	for(int rt=(mat->row);rt>=1;rt--) {
		if(*(*((mat->val)+rt-1)+(mat->col)-1)==0) 
			continue;
		for(int ct=0;ct<(mat->col)-1;ct++) 
			if(*(*((mat->val)+rt-1)+ct)!=0) 
				return 0;
	}
	return 1;
}

/*
 * sets the matrix to minimum size of one row and two columns
 */
void* mat_init() {
	struct mat tmp;
	tmp.col=2;
	tmp.row=1;
	tmp.val=malloc(sizeof(double*));
	*(tmp.val)=malloc(2*sizeof(double));
	return &tmp;
}

/*
 * destructs a matrix
 */
void mat_destr(void* mat) {
	if(!(mat==NULL)) {
		for(int rt=0;rt<(mat->row);rt++) 
			free(*((mat->val)+rt));
		free(mat->val);
		(mat->col)=0;
		(mat->row)=0;
		mat=NULL;
	}
}

/*
 * resizes the matrix to given size
 * -1 means no change here
 * warning: lastest numbers will be lost at srink
 */
void* mat_resize(void* mat,int ncol,int nrow) {
	if(!(ncol==-1 || ncol==(mat->col))) {
		double* tmp=malloc(ncol*sizeof(double));
		for(int rt=0;rt<(mat->row);rt++) {
			memcpy(tmp,*((mat->val)+rt),((mat->col)>ncol?ncol:(mat->col)));
			*((mat->val)+rt)=realloc(*((mat->val)+rt),ncol*sizeof(double));
			memcpy(*((mat->val)+rt),tmp,ncol);
		}
		(mat->col)=ncol;
		free(tmp);
	}
	if(!(nrow==-1 || nrow==(mat->row))) {
		double** tmp=malloc(nrow*sizeof(double*));
		memcpy(tmp,(mat->val),((mat->row)>nrow?nrow:(mat->row)));
		(mat->val)=realloc((mat->val),nrow*sizeof(double*));
		memcpy((mat->val),tmp,nrow);
		if(nrow>(mat->row))
			for(int rt=(mat->row);rt<nrow;rt++) 
				*((mat->val)+rt)=malloc((mat->col)*sizeof(double));
		(mat->row)=nrow;
		free(tmp);
	}
	return mat;
}

/*
 * returns value at given position
 */
double inline mat_get_val(void* mat,int r,int c) {return *(*((mat->val)+r-1)+c-1);}

/*
 * sets the value of the given position to given number
 */
void inline mat_set_val(void* mat,int r,int c,double val) {*(*((mat->val)+r-1)+c-1)=val;}

/*
 * starts solving with gauss methode
 * returns 1 if the system was solvable
 * or 0 when not
 */
int mat_solve(void* mat) {
	mat_solve_norm(mat);
	return mat_solvable(mat);
}

