#include "cgaus.h"

double** mat;
int col,row;

/*
 * returns col
 */
int inline mat_get_col() {return col;}

/*
 * returns row
 */
int inline mat_get_row() {return row;}

/*
 * multiplies every number in a line of the length col with fact
 */
void mat_mull(double* line,double fact) {
	for(int ct=0;ct<col;ct++) *(line+ct)*=fact;
}

/*
 * adds to each number in a line of the length col the specific number in the second equaly length line
 */
void mat_addl(double* tline,double* sline) {
	for(int ct=0;ct<col;ct++) *(tline+ct)+=*(sline+ct);
}

/*
 * swaps the position of the current line with a later line, 
 * that is non zero at the spos position and returns 1, 
 * if no matching line is found return 0
 */
int mat_swapl(int spos) {
	if(*(*(mat+spos-1)+spos-1)!=0) ; /*err?!*/
	double* tmpp=*(mat+spos-1);
	for(int rt=spos;rt<row;rt++) {
		if(*(*(mat+rt)+spos-1)==0)
			continue;
		*(mat+spos-1)=*(mat+rt);
		*(mat+rt)=tmpp;
		return 1;
	}
	return 0;
}

/*
 * solving the matrix to format where there is only a diagonal line of ones 
 * but the rest zeros expect last column
 */
void mat_solve_norm() {
	for(int rt=0;rt<col-1;rt++) {
		if(*(*(mat+rt)+rt)==0) 
			if(!mat_swapl(rt+1)) 
				continue;
		mat_mull(*(mat+rt),-1/(*(*(mat+rt)+rt)));
		double* tmp=malloc(col*sizeof(double));
		for(int rt2=rt+1;rt2<row;rt2++) {
			memcpy(tmp,*(mat+rt),col);
			mat_mull(tmp,*(*(mat+rt2)+rt));
			mat_addl(*(mat+rt2),tmp);
		}
		free(tmp);
	}
	for(int tr=col-1;tr>=0;tr--) {
		if(*(*(mat+tr)+tr)==0) 
			continue;
		double* tmp=malloc(col*sizeof(double));
		for(int tr2=tr-1;tr2>=0;tr2--) {
			memcpy(tmp,*(mat+tr),col);
			mat_mull(tmp,*(*(mat+tr2)+tr));
			mat_addl(*(mat+tr2),tmp);
		}
		free(tmp);
	}
}

/*
 * returns 1 if all lines with only zeros in the front ends with zero
 * else return 0
 */
int mat_solvable() {
	for(int rt=row;rt>=1;rt--) {
		if(*(*(mat+rt-1)+col-1)==0) 
			continue;
		for(int ct=0;ct<col-1;ct++) 
			if(*(*(mat+rt-1)+ct)!=0) 
				return 0;
	}
	return 1;
}

/*
 * sets the matrix to minimum size of one row and two columns
 */
void mat_init() {
	if(mat) {
		for(int r=0;r<row;r++) {
			free(*(mat+r));
		}
		free(mat);
	}
	col=2;
	row=1;
	mat=malloc(sizeof(double*));
	*mat=malloc(2*sizeof(double));
}

/*
 * extends the matrix by one column
 */
void mat_add_col() {
	double* tmp=malloc((++col)*sizeof(double));
	for(int r=0;r<row;r++) {
		memcpy(tmp,*(mat+r),col-1);
		*(mat+r)=realloc(*(mat+r),col*sizeof(double));
		memcpy(*(mat+r),tmp,col);
	}
	free(tmp);
}

/*
 * extends the matrix by one row
 */
void mat_add_row() {
	double** tmp=malloc((++row)*sizeof(double*));
	memcpy(tmp,mat,row-1);
	mat=realloc(mat,row*sizeof(double*));
	memcpy(mat,tmp,row);
	free(tmp);
	*(mat+row-1)=malloc(col*sizeof(double));
}

/*
 * returns value at given position
 */
double inline mat_get_val(int r,int c) {return *(*(mat+r-1)+c-1);}

/*
 * sets the value of the given position to given number
 */
void inline mat_set_val(int r,int c,double val) {*(*(mat+r-1)+c-1)=val;}

/*
 * starts solving with gauss methode
 * returns 1 if the system was solvable
 * or 0 when not
 */
int mat_solve() {
	mat_solve_norm();
	return mat_solvable();
}

