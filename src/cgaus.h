#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ncurses.h>

void mat_init();
int mat_solve();
void mat_add_col();
void mat_add_row();
double inline mat_get_val(int,int);
void inline mat_set_val(int,int,double);
int inline mat_get_col();
int inline mat_get_row();

