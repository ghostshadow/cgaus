#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void* mat_init(void*);
void mat_destr(void*);
int mat_solve(void*);
void* mat_resize(void*,int,int);
double mat_get_val(void*,int,int);
void mat_set_val(void*,int,int,double);
int mat_get_col(void*);
int mat_get_row(void*);

