#include "cgaus.h"

#include <stdio.h>
#include <stdlib.h>

void sif_get_mat(mxp mat) {
	int r,c;
	double tmp;
	printf("Input matrix size (<rows>x<columns>): ");
	scanf("%dx%d",&r,&c);
	mat_resize(mat,r,c);
	for(int rt=0;rt<r;rt++) 
		for(int ct=0;ct<c;ct++) {
			printf("Value of cell [%d;%d]: ",rt+1,ct+1);
			scanf("%lf",&tmp);
			mat_set_val(mat,rt,ct,tmp);
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
