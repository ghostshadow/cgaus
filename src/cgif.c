#include "cgaus.h"

#include <stdio.h>
#include <stdlib.h>
#include <curses.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

void sif_get_mat(mxp mat) {
	int r,c;
	char *tpos=NULL;
	char buffer[128];
	printf("Input matrix size (<rows>x<columns>): ");
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


int my=0,mx=0;

#define _NHELP "Help:\n"\
		"<h> to open or close help\n"\
		"<space|enter> to overwrite the actually marked cell or stop editing\n"\
		"<arrows> to move the selection\n"\
		"<tap> to calculate\n"\
		"<q> to quit\n"\
		"<n> to reset matrix to default\n"\
		"<a>[count]<l|c> to add count amount of lines|columns\n"\
		"<d>[count]<l|c> to delete count amount of lines|columns"
		/*
		 * longest line in help: 64
		 * lines in help: 9
		*/

static void nerrmsg(int y,int x,int colpair,char* errmsg) {
	int my,mx;
	getmaxyx(stdscr,my,mx);
	if(y<1 || x<1)
		return;
	int curch=0,width=0,mwidth=0,line=0;
	int* lwidth=malloc(sizeof(int));
	while(curch<(int)strlen(errmsg)) {
		if(errmsg[curch]!='\n' && (x+width+1)<mx) {
			mvaddch(y+line,x+width,errmsg[curch]);
			mwidth=(mwidth<(++width)?width:mwidth);
		}else if(errmsg[curch]=='\n' && (y+line+1)<my) {
			lwidth=realloc(lwidth,((++line)+1)*sizeof(int));
			width=0;
		}else
			return;
		lwidth[line]=width;
		curch++;
	}
	mvhline(y-1,x-1,'=',mwidth+2);
	mvhline(y+line+1,x-1,'=',mwidth+2);
	for(int l=0;l<line+1;l++) {
		mvaddch(y+l,x-1,' ');
		for(int c=lwidth[l];c<mwidth+1;c++)
			mvaddch(y+l,x+c,' ');
	}
	for(int l=0;l<line+3;l++) 
		mvchgat(y+l-1,x-1,mwidth+2,0,colpair,NULL);
	move(y+line+1,x+mwidth+1);
	refresh();
	return;
}

static void draw_headline() {
	mvprintw(0,(mx-28)/2,"Linear Equation Matrix Solver");
	mvprintw(1,(mx-17)/2,"using gaus methode");
	mvchgat(0,0,-1,A_BOLD,2,NULL);
	mvchgat(1,0,-1,0,2,NULL);
}

static void draw_helpline(char* hline1,char* hline2) {
	mvaddstr(my-2,0,hline1);
	mvaddstr(my-1,0,hline2);
	mvchgat(my-2,0,-1,0,2,NULL);
	mvchgat(my-1,0,-1,0,2,NULL);
}

static void draw_userline(char* uline) {
	move(my,0);
	clrtoeol();
	addstr(uline);
	mvchgat(my,0,-1,0,3,NULL);
}

static void draw_mat(mxp nmat,int mxy,int mxx) {
	int dsmx=mx,dsmy=my-5,dsys=2;
	mvprintw((dsmy-1)/2,(mx-17)/2,"<Matrix here ...>");
	return;
	int nmatl=mat_get_row(nmat),nmatc=mat_get_col(nmat);
	

}


static void nredraw(mxp nmat,int mxy,int mxx,char* lastcommand) {
	char* tmp2=malloc(50*sizeof(char));
	getmaxyx(stdscr,my,mx);
	my--; mx--;
	clear();
	bkgd(COLOR_PAIR(1));
	draw_mat(nmat,mxy,mxx);
	draw_headline();
	sprintf(tmp2,"Last Command: %s",lastcommand);
	draw_helpline("press h for help!",tmp2);
	draw_userline("");
	refresh();
	free(tmp2);
}

mxp* nmap;
int* mxyp,* mxxp;
char** lastcomm;

static void nresize(int sig) {
	touchwin(stdscr);
	nredraw(*nmap,*mxyp,*mxxp,*lastcomm);	
}

void stop_nif(int sig) {
	endwin();
}

int static start_nif() {
	if(!initscr())
		return 0;
	signal(SIGTERM,stop_nif);
	signal(SIGWINCH,nresize);
	cbreak();
	noecho();
	keypad(stdscr,TRUE);
	meta(stdscr,TRUE);
	if(has_colors()) {
		start_color();
		init_pair(1,COLOR_BLACK,COLOR_YELLOW);
		init_pair(2,COLOR_WHITE,COLOR_BLUE);
		init_pair(3,COLOR_WHITE,COLOR_BLACK);
		init_pair(4,COLOR_BLUE,COLOR_RED);
	}
	return 1;
}
		/*
		init_pair(1,COLOR_BLACK,COLOR_YELLOW);
		init_pair(2,COLOR_WHITE,COLOR_BLUE);
		init_pair(3,COLOR_WHITE,COLOR_BLACK);
		init_pair(4,COLOR_BLUE,COLOR_RED);
		*/

int nif() {
	if(!start_nif())
		return 0;
	getmaxyx(stdscr,my,mx);
	my--; mx--;
	mxp nmat=mat_init();
	int mxl=mat_get_row(nmat),mxc=mat_get_col(nmat),mxy=0,mxx=0;
	int nmxl=mxl,nmxc=mxc;
	int quit=0,curch=0;
	int itmp=0;
	char ctmp=0;
	double dtmp=0;
	int* tmp=malloc(50*sizeof(int));
	char* tmp2=malloc(50*sizeof(char));
	char* lastcommand=malloc(50*sizeof(char));
	nmap= &nmat;
	mxyp= &mxy;
	mxxp= &mxx;
	lastcomm= &lastcommand;
	lastcommand[0]='\0';
	tmp[0]='\0';
	tmp2[0]='\0';
	bkgd(COLOR_PAIR(1));
	nredraw(nmat,mxy,mxx,lastcommand);
	while(!quit) {
		curch=getch();
		switch(curch) {
		case 0x9:
			itmp=mat_solve(nmat);
			sprintf(lastcommand,"Solved matrix [%s]!",(itmp?"solvable":"not solvable"));
		case ' ':
		case '\n':
			sprintf(tmp2,"Choose Value for selected Cell [l%d;c%d] : ",mxy+1,mxx+1);
			draw_userline(tmp2);
			echo();
			attron(COLOR_PAIR(3));
			mvscanw(my,(int)strlen(tmp2),"%lf",&dtmp);
			attroff(COLOR_PAIR(3));
			noecho();
			mat_set_val(nmat,mxy,mxx,dtmp);
			sprintf(lastcommand,"Setting Cell Value  [l%d;c%d] to: %lf",mxy+1,mxx+1,dtmp);
			break;
		case 'Q':
		case 'q':
			quit=1;
			continue;
		case 'H':
		case 'h':
			draw_userline("Help");
			nerrmsg((my-9)/2,(mx-64)/2,2,_NHELP);
			int ch=getch(),q2=0;
			while(q2) {
				refresh();
				ch=getch();
				if(ch=='h' || ch=='q')
					q2=1;
			}
			sprintf(lastcommand,"Help");
			break;
		case 'A':
		case 'a':
			break;
			draw_userline("Adding :");
			echo();
			attron(COLOR_PAIR(3));
			itmp=1;
			mvscanw(my,8,"%s",tmp);
			attroff(COLOR_PAIR(3));
			noecho();
			if(isdigit(tmp[0]))
				sscanf((char*)tmp,"%d%c",&itmp,&ctmp);
			else
				sscanf((char*)tmp,"%c",&ctmp);
			if(ctmp=='l' || ctmp=='L') {
				nmxl=mxl+itmp;
				nmxc=mxc;
			}else if(ctmp=='c' || ctmp=='C') {
				nmxc=mxc+itmp;
				nmxl=mxl;
			}else
				break;
			mat_resize(nmat,nmxl,nmxc);
			mxl=mat_get_row(nmat); mxc=mat_get_col(nmat);
			nmxl=mxl; nmxc=mxc;
			sprintf(lastcommand,"Adding %d %s (is now: %d lines; %d columns)",itmp,(ctmp=='c' || ctmp=='C'?"columns":"lines"),mxl,mxc);
			break;
		case 'D':
		case 'd':
			break;
			draw_userline("Deleting :");
			echo();
			attron(COLOR_PAIR(3));
			itmp=1;
			mvscanw(my,8,"%s",tmp);
			attroff(COLOR_PAIR(3));
			noecho();
			if(isdigit(tmp[0]))
				sscanf((char*)tmp,"%d%c",&itmp,&ctmp);
			else
				sscanf((char*)tmp,"%c",&ctmp);
			if(ctmp=='l' || ctmp=='L') {
				nmxl=mxl-itmp;
				nmxc=mxc;
			}else if(ctmp=='c' || ctmp=='C') {
				nmxc=mxc-itmp;
				nmxl=mxl;
			}else
				break;
			mat_resize(nmat,nmxl,nmxc);
			mxl=mat_get_row(nmat); mxc=mat_get_col(nmat);
			nmxl=mxl; nmxc=mxc;
			sprintf(lastcommand,"Deleting %d %s (is now: %d lines; %d columns)",itmp,(ctmp=='c' || ctmp=='C'?"columns":"lines"),mxl,mxc);
			break;
		case KEY_UP:
			if(mxy>0) {
				mxy--;	
				sprintf(lastcommand,"Moved to [l%d;c%d] with actual value: %lf",mxy+1,mxx+1,mat_get_val(nmat,mxy,mxx));
			}
			break;
		case KEY_DOWN:
			if(mxy<mxl-1) {
				mxy++;
				sprintf(lastcommand,"Moved to [l%d;c%d] with actual value: %lf",mxy+1,mxx+1,mat_get_val(nmat,mxy,mxx));
			}
			break;
		case KEY_LEFT:
			if(mxx>0) {
				mxx--;
				sprintf(lastcommand,"Moved to [l%d;c%d] with actual value: %lf",mxy+1,mxx+1,mat_get_val(nmat,mxy,mxx));
			}
			break;
		case KEY_RIGHT:
			if(mxx<mxc-1) {
				mxx++;
				sprintf(lastcommand,"Moved to [l%d;c%d] with actual value: %lf",mxy+1,mxx+1,mat_get_val(nmat,mxy,mxx));
			}
			break;
		case 'P':
		case 'p':
			sprintf(lastcommand,"Size: [l%d;c%d]; Position: [l%d;c%d]; Cell Value: %lf",mxl,mxc,mxy+1,mxx+1,mat_get_val(nmat,mxy,mxx));
			break;
		}
		nredraw(nmat,mxy,mxx,lastcommand);
	}
	stop_nif(0);
	mat_destr(nmat);
	signal(SIGTERM,SIG_DFL);
	signal(SIGWINCH,SIG_DFL);
	return 1;
}





