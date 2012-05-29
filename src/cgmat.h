typedef struct mat* mxp;

mxp mat_init();
void mat_destr(mxp);
int mat_solve(mxp);
mxp mat_resize(mxp,int,int);
double mat_get_val(mxp,int,int);
void mat_set_val(mxp,int,int,double);
int mat_get_col(mxp);
int mat_get_row(mxp);

