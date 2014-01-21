#include "cgaus.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void sif_get_mat(mxp mat) {
	int r,c;
	char *tpos=NULL;
	char buffer[128];
	printf("GAUSS LINEAR EQUATION SOLVER\nInput matrix size including result vector as additional column \n(<rows>x<columns>): ");
	fgets(&buffer[0],128,stdin);
	r=strtol(&buffer[0],&tpos,10);
	c=strtol(++tpos,&tpos,10);
	mat_resize(mat,r,c);
	for(int rt=0;rt<r;rt++) 
		for(int ct=0;ct<c;ct++) {
			memset(&buffer[0],0,128);
			printf("Value of cell [%d;%d]: ",rt+1,ct+1);
			fgets(&buffer[0],128,stdin);
			mat_set_val(mat,rt,ct,strtod(&buffer[0],NULL));
		}
	putchar('\n');
}

void sif_print_mat(mxp mat) {
	for(int rt=0;rt<mat_get_row(mat);rt++) {
		for(int ct=0;ct<mat_get_col(mat);ct++) {
			printf("%lf ",mat_get_val(mat,rt,ct));
		}
		putchar('\n');
	}
	putchar('\n');
}

