#include <stdio.h>
#include <math.h>
#include <time.h>

#define NX 500											/****************************************************/
#define NY 50											/*...has to be consistent with the Global NS mesh...*/
#define NZ 100											/****************************************************/

int i,j,k;
double L,B,H,dx,dy,dz,px[NX+10],py[NY+10],pz[NZ+10],cx[NX+10],cy[NY+10],cz[NZ+10];
double lx[NX+10],ly[NY+10],lz[NZ+10],rx,ry,rz,exsum;
int IMG=4;											//number of imaginary cells - has to be consistent with "3-D_TWO-PHASE-NSE.CPP"
double f[NX+10][NY+10][NZ+10]={0.0};

FILE *f1;

double pi=3.141592653589793238462643383279502884197;
double a_cyl=10.0,b_cyl=2.0,r=0.5;								//cylinder centered at (x=a_cyl,y=b_cyl) having a radius of "r"

main()
{
clock_t tt;
tt=clock();

printf("This is a code to PATCH the FFD VOF-field for a circular cylinder\n\n");

L=20.0;												//length of the NS domain
B=4.0;												//width of the NS domain
H=2.0;												//height of the NS domain

rx=1.0;												//mesh-stretching parameter along x
ry=1.0;												//mesh-stretching parameter along y
rz=1.0;												//mesh-stretching parameter along z

//computing cell-sizes along the x-direction
exsum=1.0;
for (i=1;i<=(NX-1);i++)
    exsum=exsum+pow(rx,i);

lx[IMG+1]=L/exsum;										//first length in x
for (i=IMG+2;i<=NX+IMG;i++)
    lx[i]=lx[i-1]*rx;

//computing cell-sizes along the y-direction
exsum=1.0;
for (j=1;j<=(NY-1);j++)
    exsum=exsum+pow(ry,j);

ly[IMG+1]=B/exsum;										//first length in y
for (j=IMG+2;j<=NY+IMG;j++)
    ly[j]=ly[j-1]*ry;

//computing cell-sizes along the z-direction
exsum=1.0;
for (k=1;k<=(NZ-1);k++)
    exsum=exsum+pow(rz,k);

lz[IMG+1]=H/exsum;										//first length in z
for (k=IMG+2;k<=NZ+IMG;k++)
    lz[k]=lz[k-1]*rz;

//calculating the coordinates of cell faces
px[IMG+1]=0.0;
py[IMG+1]=0.0;
pz[IMG+1]=0.0;
for (i=IMG+2;i<=NX+IMG+1;i++)
    px[i]=px[i-1]+lx[i-1];
for (j=IMG+2;j<=NY+IMG+1;j++)
    py[j]=py[j-1]+ly[j-1];
for (k=IMG+2;k<=NZ+IMG+1;k++)
    pz[k]=pz[k-1]+lz[k-1];

//calculating the coordinates of cell centers
cx[IMG+1]=0.5*lx[IMG+1];
cy[IMG+1]=0.5*ly[IMG+1];
cz[IMG+1]=0.5*lz[IMG+1];
for (i=IMG+2;i<=NX+IMG+1;i++)
    cx[i]=cx[i-1]+0.5*(lx[i-1]+lx[i]);
for (j=IMG+2;j<=NY+IMG+1;j++)
    cy[j]=cy[j-1]+0.5*(ly[j-1]+ly[j]);
for (k=IMG+2;k<=NZ+IMG+1;k++)
    cz[k]=cz[k-1]+0.5*(lz[k-1]+lz[k]);

//pixelization routine for patching
for (i=IMG+1;i<=NX+IMG;i++)
{
    for (j=IMG+1;j<=NY+IMG;j++)
    {
        for (k=IMG+1;k<=NZ+IMG;k++)								//cylinder spans the entire height of the domain
        {
            /*...Checking whether ANY of the 8 hexahedron vertices fall INSIDE the cylinder...*/
	    if ((pow((px[i]-a_cyl),2.0)+pow((py[j]-b_cyl),2.0)<r*r) || (pow((px[i+1]-a_cyl),2.0)+pow((py[j]-b_cyl),2.0)<r*r) || (pow((px[i]-a_cyl),2.0)+pow((py[j+1]-b_cyl),2.0)<r*r) || (pow((px[i+1]-a_cyl),2.0)+pow((py[j+1]-b_cyl),2.0)<r*r))
                {
                    /*...Checking whether ALL of the 8 hexahedron vertices fall INSIDE the cylinder...*/
		    if ((pow((px[i]-a_cyl),2.0)+pow((py[j]-b_cyl),2.0)<r*r) && (pow((px[i+1]-a_cyl),2.0)+pow((py[j]-b_cyl),2.0)<r*r) && (pow((px[i]-a_cyl),2.0)+pow((py[j+1]-b_cyl),2.0)<r*r) && (pow((px[i+1]-a_cyl),2.0)+pow((py[j+1]-b_cyl),2.0)<r*r))
                        f[i][j][k]=1.0;								//pure cell (VOF_SOLID=1)
                    else									//mixed cell (0<VOF_SOLID<1)
                    {
                        int l,m,n,VOFcount=0,SUBcells=200;
                        double dxx,dyy,dzz;
                        dxx=lx[i]/SUBcells;
                        dyy=ly[j]/SUBcells;
                        dzz=lz[k]/SUBcells;
			//generating a SUBcells^3 sub-cell mesh (face-centered)
			double pxx[SUBcells+10],pyy[SUBcells+10],pzz[SUBcells+10];
			    pxx[1]=px[i];
			    pyy[1]=py[j];
			    pzz[1]=pz[k];
			for (l=2;l<=SUBcells+1;l++)
			    pxx[l]=pxx[l-1]+dxx;
			for (m=2;m<=SUBcells+1;m++)
			    pyy[m]=pyy[m-1]+dyy;
			for (n=2;n<=SUBcells+1;n++)
			    pzz[n]=pzz[n-1]+dzz;

                            for (l=1;l<=SUBcells;l++)
                            {
                                for (m=1;m<=SUBcells;m++)
                                {
                                    for (n=1;n<=SUBcells;n++)
                                    {
				    /*...Checking whether ALL of the 8 SUBcell hexahedron vertices fall INSIDE the cylinder...*/
				    if ((pow((pxx[l]-a_cyl),2.0)+pow((pyy[m]-b_cyl),2.0)<r*r) && (pow((pxx[l+1]-a_cyl),2.0)+pow((pyy[m]-b_cyl),2.0)<r*r) && (pow((pxx[l]-a_cyl),2.0)+pow((pyy[m+1]-b_cyl),2.0)<r*r) && (pow((pxx[l+1]-a_cyl),2.0)+pow((pyy[m+1]-b_cyl),2.0)<r*r))
                                        VOFcount=VOFcount+1;
                                    }
                                }
                            }
                        f[i][j][k]=(double)VOFcount/(double)(pow(SUBcells,3.0));
                        printf("The code is here: i=%d\tj=%d\tk=%d\tf=%f\tVOFcount=%d\n",i,j,k,f[i][j][k],VOFcount);
                    }
                }
            else										//pure cell (VOF_SOLID=0)
                continue;
        }
    }
}

//writing the initialized/patched VOF field to file - TO BE READ BY IITM-RANS3D
f1=fopen("CYLINDER.dat","w");
for (i=IMG+1;i<=NX+IMG;i++)
    for (j=IMG+1;j<=NY+IMG;j++)
        for (k=IMG+1;k<=NZ+IMG;k++)
        fprintf(f1,"%f\t%f\t%f\t%0.10f\n",cx[i],cy[j],cz[k],f[i][j][k]);
fclose(f1);

//writing the initialized/patched VOF field to file - FOR CHECKING QUALITY OF THE PATCH IN TECPLOT/PARAVIEW
f1=fopen("VOFINIT_tecplot.dat","w");
fprintf (f1,"\nVARIABLES = ""X"", ""Y"", ""Z"", ""VOF""\n");
fprintf (f1,"ZONE I=%d, J=%d, K=%d, ZONETYPE=ORDERED, DATAPACKING=POINT\n",NZ,NY,NX);
for (i=IMG+1;i<=NX+IMG;i++)
    for (j=IMG+1;j<=NY+IMG;j++)
        for (k=IMG+1;k<=NZ+IMG;k++)
        fprintf(f1,"%f\t%f\t%f\t%0.10f\n",cx[i],cy[j],cz[k],f[i][j][k]);
fclose(f1);

tt=clock()-tt;
printf("The computation time is %f seconds\n",(float)tt/CLOCKS_PER_SEC);
}
