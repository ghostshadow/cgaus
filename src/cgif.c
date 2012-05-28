#include "cgaus.h"

#include <stdio.h>
#include <stdlib.h>

void sif_get_mat(void* mat) {
	int c,r;
	double tmp;
	printf("Input matrix size (<columns>x<rows>): ");
	scanf("%dx%d",&c,&r);
	mat_resize(mat,c,r);
	for(int rt=0;rt<r;rt++) 
		for(int ct=0;ct<c;ct++) {
			printf("Value of cell [%d;%d]: ",ct+1,rt+1);
			scanf("%lf",&tmp);
			mat_set_val(mat,ct,rt,tmp);
		}
	putchar('\n');
}

void sif_print_mat(void* mat) {
	for(int rt=0;rt<mat_get_row(mat);rt++) {
		for(int ct=0;ct<mat_get_col(mat);ct++) {
			printf("%lf ",mat_get_val(mat,ct,rt));
		}
		putchar('\n');
	}
	putchar('\n');
}
