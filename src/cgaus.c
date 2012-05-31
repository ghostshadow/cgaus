#include "cgaus.h"

#include <stdio.h>

void sif() {
	mxp mat=mat_init();
	sif_get_mat(mat);
	sif_print_mat(mat);
	printf("\nMatrix is [%s]:\n\n",(mat_solve(mat)?"solvable":"not solvable"));
	sif_print_mat(mat);
	mat_destr(mat);
}

int main(int argc, char** argv) {
	sif();
	//init
	//mainloop
	//bla
	return 0;
}


