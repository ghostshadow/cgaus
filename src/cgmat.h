#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void* mat_init();
void mat_destr(void*);
int mat_solve(void*);
void* mat_resize(void*,int,int);
double inline mat_get_val(void*,int,int);
void inline mat_set_val(void*,int,int,double);
int inline mat_get_col(void*);
int inline mat_get_row(void*);

