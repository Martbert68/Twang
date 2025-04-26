#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#include "martin.h"
#define BRIDGE 50
#define DOWN 250
#define UP 250

/* here are our X variables */
Display *dis;

int screen,len;
double scale;
double mass[6];
int fret[6];
unsigned char *x_buffer;
Window win;
GC gc;
XImage *x_image;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();
int f_col,f_row,f_off,dot_x,dot_y;
static double pp; 


struct node { double x;  double y; double z; double a; double b; double c; double xv; double yv; double zv; double av; double bv; double cv;};


void disp (unsigned char *,int,int);

void usage ()
{
	printf("usage: font filename threshold [20-40 ish] star,framestcode [65A 97a]\n");
	exit (1);
}

int qlott(unsigned char *image, int xp, int yp, int r, int g, int b, int xs, int ys)
{
	int x,y;
	if (xs<1){ xs=1;}
	if (ys<1){ ys=1;}
	for (x=0;x<xs;x++)
	{
		for (y=0;y<ys;y++)
		{
			plott(image,xp+x,yp+y,r,g,b);
		}
	}
}	

int load_font(char *name,unsigned int *font_n)
{

  FILE * font_file;             /* source file */
  char *line,*test;
  size_t len = 0;
  ssize_t read;
  int count;
  count=0;

  line=(char *)malloc(100);


  if ((font_file = fopen(name, "r")) == NULL) {
    fprintf(stderr, "can't open %s\n", name);
    return 0;}


    int bit_map,no_char,scanning;
    no_char=1;
    bit_map=0;
    scanning=1;

    while ((read = getline(&line, &len, font_file)) != -1) {
        //printf("Retrieved line of length %zu:\n", read);
        //printf("%s", line);
        if (no_char)
        {
                if (!strncmp(line,"STARTCHAR",9))
                {
                        no_char=0;
                        //printf("start %s \n",line);
                }
        }else{
                if (!strncmp(line,"ENDCHAR",7))
                {
                        no_char=1;
                        //printf("End %s \n",line);
                        bit_map=0;
                }
                if (bit_map)
                {
                        unsigned int bit;
                        sscanf(line, "%X", &bit);
                        //printf("Bit %s %d\n",line,bit);
                        //DisplayBinary(bit);
                        //printf("\n");
                        font_n[count]=bit;
                        count++;
                }
                if (!strncmp(line,"BITMAP",6))
                {
                        bit_map=1;
                }
                if (!strncmp(line,"BBX",3) && scanning)
                {
                        int r,c,la,off;
                        char bbx[10];
                        sscanf(line,"%s %d %d %d %d",bbx,&c,&r,&la,&off);
                        printf ("c %d r %d off %d\n",c,r,-off);
                        if (off==-2){off=0;}
                        if (c>1 && c<50 && r>1 && r<50){
                                scanning=0;
                                f_col=c;f_row=r;f_off=-off+2;
                                printf ("off %d\n",f_off);
                        }
	          }
                }
        }

  fclose(font_file);
  free (line);

 if (!scanning)
  {
        return count;
  }else{
          return 0;
  }

}

int draw_letter(unsigned char *image1, int x,int y, int red, int green, int blue, unsigned int *font_n, int letter_no, float theta)
{
        int r,c,xw,yw,i;

        float dx,dy;

        dx=cos(theta);
        dy=sin(theta);

        if (letter_no<0){ return 0;}

        for (r=0;r<f_row;r++)
        {
                unsigned int row_v;
                unsigned int point;
                point=1<<(f_col+f_off);
                row_v=font_n[(f_row*letter_no)+(r)];
                for (c=0;c<f_col;c++)
                {
                        if ( row_v & point )
                        {
                                int xp,yp;
                                for (xw=0;xw<dot_x;xw++)
                                {
                                        for (yw=0;yw<dot_y;yw++)
                                        {
                                                xp=(x+xw+((float)c*dx*dot_x)+((float)r*dy*dot_y));
                                                if (xp>=X_SIZE){continue;}
                                                if (xp<0){continue;}
                                                yp=(y+yw-((float)c*dy*dot_y)+((float)r*dx*dot_x));
                                                if (yp>=Y_SIZE){continue;}
                                                if (yp<0){continue;}
                                                image1[(xp*3)+(yp*3*X_SIZE)]=red;
                                                image1[1+(xp*3)+(yp*3*X_SIZE)]=green;
                                                image1[2+(xp*3)+(yp*3*X_SIZE)]=blue;

                                                //draw_point (zframe, xp, yp, red,green,blue,power);
                                        }
                                }
                        }
                        point=point >>1;
                        //printf ("%d point %d row_v \n",point,row_v);
                }
        }

        return 1;
}


int draw_word (char *string, int *xp, int *yp, unsigned int *font_n, float ang, float wang,int r, int g, int b, unsigned char *image5)
{
        int p;
        p=0;
        char c;
        c=string[p];
        float xd,yd;
        xd=sin(wang);
        yd=cos(wang);
        while (c>0)
        {
                printf ("%c",c);
                draw_letter(image5,*xp+(p*dot_x*yd*f_col),*yp-(p*dot_y*xd*f_row),r,g,b,font_n,c-31,ang);
                p++;
                c=string[p];
        }
        *xp=*xp+(p*dot_x*yd*f_row);
        *yp=*yp-(p*dot_y*xd*f_col);
        printf("\n");
}

void init_n(struct node **n ,struct node **m)
{
	int i;
       	*(n)=(struct node *)malloc(sizeof (struct node)); //nodes 
       	*(m)=(struct node *)malloc(sizeof (struct node)); //nodes 
       	*(n+len+1)=(struct node *)malloc(sizeof (struct node)); //nodes 
       	*(m+len+1)=(struct node *)malloc(sizeof (struct node)); //nodes 
	for (i=0;i<len+2;i++)
	{
        	*(n+i)=(struct node *)malloc(sizeof (struct node)); //nodes 
        	*(m+i)=(struct node *)malloc(sizeof (struct node)); //nodes 
		//n[i]->x=0; n[i]->y=0; n[i]->z=0; n[i]->xv=0; n[i]->yv=0; n[i]->m=2+(4*(double)i/(double)len);
		n[i]->x=0; n[i]->y=0; n[i]->z=0; n[i]->xv=0; n[i]->yv=0; 
		int s;
	}

	double mm;
	mm=25;
	mass[0]=mm;
	mass[1]=mm/pow(pow(pp,5),2);
	mass[2]=mm/pow(pow(pp,10),2);
	mass[3]=mm/pow(pow(pp,15),2);
	mass[4]=mm/pow(pow(pp,19),2);
	mass[5]=mm/pow(pow(pp,24),2);

}

void do_n(struct node **n ,struct node **m)
{
	int i;
	double xa,ya,za,aa,ba,ca,d,bb,cc,xx,yy,zz;

	d=0.9999965;
	d=0.99999;
	//d=0.9999997;
	//d=0.9999999;

	for (i=1;i<len+1;i++)
	{
		// a=F/m

		/*if (i==BRIDGE)
		{
		xa=((n[i-1]->x+n[i+1]->x)-(2*n[i]->x))/30;
		ya=((n[i-1]->y+n[i+1]->y)-(2*n[i]->y))/30;
		za=((n[i-1]->z+n[i+1]->z)-(2*n[i]->z))/30;
		aa=((n[i-1]->a+n[i+1]->a)-(2*n[i]->a))/30;
		ba=((n[i-1]->b+n[i+1]->b)-(2*n[i]->b))/30;
		ca=((n[i-1]->c+n[i+1]->c)-(2*n[i]->c))/30;
		} else{ */
		xa=((n[i-1]->x+n[i+1]->x)-(2*n[i]->x))/mass[0];
		ya=((n[i-1]->y+n[i+1]->y)-(2*n[i]->y))/mass[1];
		za=((n[i-1]->z+n[i+1]->z)-(2*n[i]->z))/mass[2];
		aa=((n[i-1]->a+n[i+1]->a)-(2*n[i]->a))/mass[3];
		ba=((n[i-1]->b+n[i+1]->b)-(2*n[i]->b))/mass[4];
		ca=((n[i-1]->c+n[i+1]->c)-(2*n[i]->c))/mass[5];
		//}
		
		m[i]->yv=(n[i]->yv*d)+ya;
		m[i]->xv=(n[i]->xv*d)+xa;
		m[i]->zv=(n[i]->zv*d)+za;
		m[i]->av=(n[i]->av*d)+aa;
		m[i]->bv=(n[i]->bv*d)+ba;
		m[i]->cv=(n[i]->cv*d)+ca;

		m[i]->x=n[i]->x+m[i]->xv; 
		m[i]->y=n[i]->y+m[i]->yv; 
		m[i]->z=n[i]->z+m[i]->zv; 
		m[i]->a=n[i]->a+m[i]->av; 
		m[i]->b=n[i]->b+m[i]->bv; 
		m[i]->c=n[i]->c+m[i]->cv; 


	}
	for (i=1;i<len+1;i++)
	{
		/*if (i==BRIDGE)
		{
		xa=((m[i-1]->x+m[i+1]->x)-(2*m[i]->x))/30;
		ya=((m[i-1]->y+m[i+1]->y)-(2*m[i]->y))/30;
		za=((m[i-1]->z+m[i+1]->z)-(2*m[i]->z))/30;
		aa=((m[i-1]->a+m[i+1]->a)-(2*m[i]->a))/30;
		ba=((m[i-1]->b+m[i+1]->b)-(2*m[i]->b))/30;
		ca=((m[i-1]->c+m[i+1]->c)-(2*m[i]->c))/30;
		}else{ */
		xa=((m[i-1]->x+m[i+1]->x)-(2*m[i]->x))/mass[0];
		ya=((m[i-1]->y+m[i+1]->y)-(2*m[i]->y))/mass[1];
		za=((m[i-1]->z+m[i+1]->z)-(2*m[i]->z))/mass[2];
		aa=((m[i-1]->a+m[i+1]->a)-(2*m[i]->a))/mass[3];
		ba=((m[i-1]->b+m[i+1]->b)-(2*m[i]->b))/mass[4];
		ca=((m[i-1]->c+m[i+1]->c)-(2*m[i]->c))/mass[5];
		//}

		n[i]->xv=m[i]->xv+xa;
		n[i]->yv=m[i]->yv+ya;
		n[i]->zv=m[i]->zv+za;
		n[i]->av=m[i]->av+aa;
		n[i]->bv=m[i]->bv+ba;
		n[i]->cv=m[i]->cv+ca;

		n[i]->x=m[i]->x+n[i]->xv; 
		n[i]->y=m[i]->y+n[i]->yv; 
		n[i]->z=m[i]->z+n[i]->zv; 
		n[i]->a=m[i]->a+n[i]->av; 
		n[i]->b=m[i]->b+n[i]->bv; 
		n[i]->c=m[i]->c+n[i]->cv; 

		if (i==fret[0]){ n[i]->x=0;n[i]->xv=0;}
		if (i==fret[1]){ n[i]->y=0;n[i]->yv=0;}
		if (i==fret[2]){ n[i]->z=0;n[i]->zv=0;}
		if (i==fret[3]){ n[i]->a=0;n[i]->av=0;}
		if (i==fret[4]){ n[i]->b=0;n[i]->bv=0;}
		if (i==fret[5]){ n[i]->c=0;n[i]->cv=0;}
	}



	aa=n[BRIDGE]->av;
	bb=n[BRIDGE]->bv;
	cc=n[BRIDGE]->cv;
	xx=n[BRIDGE]->xv;
	yy=n[BRIDGE]->yv;
	zz=n[BRIDGE]->zv;


	n[BRIDGE]->xv+=(-xx+yy+zz+aa+bb+cc)/1200;
	n[BRIDGE]->yv+=(xx-yy+zz+aa+bb+cc)/1200;
	n[BRIDGE]->zv+=(xx+yy-zz+aa+bb+cc)/1200;
	n[BRIDGE]->av+=(xx+yy+zz-aa+bb+cc)/1200;
	n[BRIDGE]->bv+=(xx+yy+zz+aa-bb+cc)/1200;
	n[BRIDGE]->cv+=(xx+yy+zz+aa+bb-cc)/1200;

}


void draw_n (unsigned char *image, struct node **n,int frame)
{
	int i;
	for (i=0;i<len+2;i++)
	{

	//qlott(image, n[i]->x+((X_SIZE)/2), n[i]->y+(Y_SIZE/2), 255*(sin(i)), 255*(sin(1.3*(double)i)), 255*(cos((double)i*3.5)), n[i]->m,n[i]->m);
	//
	//
	int r,g,b;

	r=128+(127*(sin((double)(i-(3*frame))/10))); g=128+(127*(cos((double)(i+(4*frame))/37))); b=128+(127*(sin((double)(i+frame)/59)));
	qlott(image, i*X_SIZE/len, (n[i]->x/4)+(Y_SIZE/7), r,g,b, 1,12);
	r=128+(127*(sin((double)(i-(3*frame))/11))); g=128+(127*(cos((double)(i+(4*frame))/39))); b=128+(127*(sin((double)(i+frame)/51)));
	qlott(image, i*X_SIZE/len, (n[i]->y/4)+(2*Y_SIZE/7), r,g,b, 1,10);
	r=128+(127*(sin((double)(i-(3*frame))/19))); g=128+(127*(cos((double)(i+(4*frame))/31))); b=128+(127*(sin((double)(i+frame)/48)));
	qlott(image, i*X_SIZE/len, (n[i]->z/4)+(3*Y_SIZE/7), r,g,b, 1,8);
	r=128+(127*(sin((double)(i-(3*frame))/12))); g=128+(127*(cos((double)(i+(4*frame))/32))); b=128+(127*(sin((double)(i+frame)/41)));
	qlott(image, i*X_SIZE/len, (n[i]->a/4)+(4*Y_SIZE/7), r,g,b, 1,6);
	r=128+(127*(sin((double)(i-(3*frame))/13))); g=128+(127*(cos((double)(i+(4*frame))/38))); b=128+(127*(sin((double)(i+frame)/32)));
	qlott(image, i*X_SIZE/len, (n[i]->b/4)+(5*Y_SIZE/7), r,g,b, 1,4);
	r=128+(127*(sin((double)(i-(3*frame))/39))); g=128+(127*(cos((double)(i+(4*frame))/33))); b=128+(127*(sin((double)(i+frame)/21)));
	qlott(image, i*X_SIZE/len, (n[i]->c/4)+(6*Y_SIZE/7), r,g,b, 1,2);
	}
	for (i=0;i<6;i++)
	{
		qlott(image, fret[i]*X_SIZE/len, (i+1)*Y_SIZE/7, 255,255,255, 30,30);
	}
}

int main(int argc,char *argv[])
{
	unsigned char *image2,*image3,*image4,*image5,*image6;
	int i,loop,frame;
	FILE *list;
        unsigned int *font_n;
	struct node **n,**m;
	short *wav;

	dot_x=4;dot_y=4;
	pp=pow(2,(double)1/(double)12);

	len=1141;


	for  (i=0;i<6;i++){fret[i]=len;}

	n = (struct node**)malloc(len * sizeof(struct node*)); 
	m = (struct node**)malloc(len * sizeof(struct node*)); 

        font_n=(unsigned int *)malloc(sizeof (unsigned int)*9*5000*15); // font
	init_n(n,m);



								    
	//load_font("9x15.bdf",font_n);

	char junk[30];



	init_x();

	frame=0;
	int rate,dur,total,frames;
	int pluck[1000][7];


	int a,b;
	for (a=0;a<1000;a++)
	{
	for (b=0;b<7;b++)
	{
		pluck[a][b]=-1;
	}
	}



  if ((list= fopen("tune.lst", "r")) == NULL) {
    fprintf(stderr, "can't open %s\n", "tune.lst");
    return 0;}


	int nc,pc,ll;
	char line[300];
	pc=0;
    while (fgets(line,300,list) != NULL) {
	    int k;
	    for (k=0;k<7;k++)
	    {
		    if (k>=strlen(line)-1){ pluck[pc][k]=-1; continue;}
		    if (line[k]==' '){ pluck[pc][k]=-1; continue;}
		    if (line[k]>47 && line[k]<58){ pluck[pc][k]=(int)line[k]-48;}
		    else if(line[k]>64 && line[k]<88){ pluck[pc][k]=(int)line[k]-55;}
		    else{printf ("EE %c EE\n",line[k]);exit(1);}
	    }
	    printf ("%s",line);
	    printf ("%d %d %d %d %d %d   %d\n",pluck[pc][0],pluck[pc][1],pluck[pc][2],pluck[pc][3],pluck[pc][4],pluck[pc][5],pluck[pc][6]);
	pc++;
    }

        printf("I found %d notes\n", pc);
	nc=pc;

	rate=48000;
	dur=300;

	frames=30;
	

        wav=(short *)malloc(sizeof (short)*2*rate*dur); // font


	total=0;
	frame=0;

	pc=0;
	int ss;
	ss=0;

        image2=(unsigned char *)malloc(sizeof (char)*X_SIZE*Y_SIZE*3); // disp buffer
								       
	int go;
	go=1;
	while (go)
	{
	int j;

	if (total%(2*rate)==0){ss=400;} 
	if (total%(rate/4)==0 || (pluck[pc][6]==1 && total%(rate/8)==0))
	{
	
	       	//for (j=0;j<6;j++){fret[j]=len;}	
	 	if (pluck[pc][5]>-1){ n[DOWN]->xv=(400+ss)/sqrt(mass[0]); fret[0]=len/pow(pp,pluck[pc][5]);}
	 	if (pluck[pc][4]>-1){ n[DOWN]->yv=(400+ss)/sqrt(mass[1]); fret[1]=len/pow(pp,pluck[pc][4]);}
	 	if (pluck[pc][3]>-1){ n[DOWN]->zv=(400+ss)/sqrt(mass[2]); fret[2]=len/pow(pp,pluck[pc][3]);}
	 	if (pluck[pc][2]>-1){ n[DOWN]->av=(400+ss)/sqrt(mass[3]); fret[3]=len/pow(pp,pluck[pc][2]);}
	 	if (pluck[pc][1]>-1){ n[DOWN]->bv=(400+ss)/sqrt(mass[4]); fret[4]=len/pow(pp,pluck[pc][1]);}
	 	if (pluck[pc][0]>-1){ n[DOWN]->cv=(400+ss)/sqrt(mass[5]); fret[5]=len/pow(pp,pluck[pc][0]);}
		printf("Note %d\n",pc);
		pc++;
		ss=0;
		if (pc>nc){go=0;}
	}

	for (j=0;j<10;j++)
	{
	do_n(n,m);
	}
	if (total%(rate/frames)==0)
	{
	blurt(image2,1);
	//clear(image2,0);
	draw_n(image2,n,frame);
	disp(image2,frame,1);
	frame++;
	}
	double s,t;
	double   ja,jp;
	s=0;t=0;
	ja=0;jp=1;s+=n[BRIDGE]->x*ja; t+=n[BRIDGE]->x*jp;
	ja=0.8;jp=0.2;s+=n[BRIDGE]->y*ja; t+=n[BRIDGE]->y*jp;
	ja=0.6;jp=0.4;s+=n[BRIDGE]->z*ja; t+=n[BRIDGE]->z*jp;
	ja=0.4;jp=0.6;s+=n[BRIDGE]->a*ja; t+=n[BRIDGE]->a*jp;
	ja=0.2;jp=0.8;s+=n[BRIDGE]->b*ja; t+=n[BRIDGE]->b*jp;
	ja=1;jp=0;s+=n[BRIDGE]->c*ja; t+=n[BRIDGE]->c*jp;
	s*=127680/Y_SIZE; t*=127680/Y_SIZE;
	if (s>32700){s=32700;} if (s<-32700){s=-32700;}
	if (t>32700){t=32700;} if (t<-32700){t=-32700;}
	wav[2*total]=s;
	wav[(2*total)+1]=t;
	total++;
	}

	save_wav(wav,"twang.wav",rate, 2, total*2 );


	scanf("%c",junk);


	close_x();

	exit(0);
}	

void disp (unsigned char *image2,int fram,int ab)
{
	int x,y;
	char *input;
	input=malloc(300);


       	for (y=0;y<Y_SIZE;y++)
       	{
               	int p=y*X_SIZE*3;
               	int XYP=X_SIZE*4*y;
               	for (x=0;x<X_SIZE;x++)
               	{
			int xpoint;
			int X_POINT;
			X_POINT=XYP+(4*x);
			xpoint=(x*3)+(p);

			x_buffer[X_POINT+2]=image2[xpoint];
			x_buffer[X_POINT+1]=image2[xpoint+1];
			x_buffer[X_POINT]=image2[xpoint+2];
                }
        }
	XPutImage(dis, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	sprintf(input,"./jpegs/pl%05d.jpg",fram);
	if (ab){jayit(image2,X_SIZE, Y_SIZE, input);}
	free (input);
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
        unsigned long black,white;

        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //y_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //z_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        dis=XOpenDisplay((char *)0);
        screen=DefaultScreen(dis);
        black=BlackPixel(dis,screen),
        white=WhitePixel(dis,screen);
        win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        XSetStandardProperties(dis,win,"image","images",None,NULL,0,NULL);
        gc=XCreateGC(dis, win, 0,0);
        XSetBackground(dis,gc,black); XSetForeground(dis,gc,white);
        XClearWindow(dis, win);
        XMapRaised(dis, win);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(dis, 0);
        x_image=XCreateImage(dis, visual, DefaultDepth(dis,DefaultScreen(dis)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(dis, gc);
        XDestroyWindow(dis,win);
        XCloseDisplay(dis);
        exit(1);
};

void redraw() {
        XClearWindow(dis, win);
};

