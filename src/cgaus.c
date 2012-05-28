#include "cgaus.h"

void sif() {
	void* mat=mat_init();
	sif_get_mat(mat);
	mat_solve(mat);
	sif_print_mat(mat);
}

int main(int argc, char** argv) {
	sif();
	//init
	//mainloop
	//bla
	return 0;
}


