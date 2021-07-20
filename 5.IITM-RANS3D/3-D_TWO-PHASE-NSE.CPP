//...THIS IS STAGE 3 OF THE DEVELOPMENT OF THE 3D HYBRID FNPT-NS COUPLED CODE...//
//...INSTANTANEOUS NSE DISCRETIZED IN THE FVM FRAMEWORK ON A STAGGERED GRID...//
//...OPERATOR-SPLIT CICSAM COUPLED WITH PARKER AND YOUNGS METHOD-BASED INTERFACE ANGLE CALCULATION...//
//...ALTERNATE SWEEPING AND CONSERVATIVE REDISTRIBUTION ALGORITHM FOR VOLUME PRESERVATION...//
//...INVERSE WEIGHTED DISTANCE-BASED HARMONIC AVERAGING FOR CALCULATION OF VISCOSITY...//

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "functions.h"

////////////////////////////////////////////////////FUNCTIONS END\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\...................
/////////////////////////////////////////////////////MAIN BEGINS\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\...................
//////////////////////////////////////////THE SOLVER LOOP IS WITHIN MAIN()\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\...................
main()
{
clock_t tt;
tt=clock();

/*...FLUID PROPERTIES...*/
rho_liquid=1000.0;										//density of water
mu_liquid=0.001137;										//dynamic viscosity of water
rho_gas=1.226;											//density of air
mu_gas=0.0000178;										//dynamic viscosity of air
rho_solid=1000.0;										//density of the "SOLID" phase
mu_solid=1.137;											//viscosity of the "SOLID" phase

/*...ACCELERATION DUE TO GRAVITY...*/
g=9.81;

/*...NUMBER OF IMAGINARY CELLS...*/
IMG=4;												//do not change - has to be consistent with "VOF_INIT.c"

/*...Relaxation parameter and error limit for GSSOR...*/
ap=0.0;												//do not change
EOPC_ERROR=0.0000001;										//RMS convergence criterion for pressure-correction equation

/*...dimensions of the domain...*/
L=20.0;												//length of the domain
B=4.0;												//width of the domain
H=2.0;												//height of the domain

/*...still water depth...*/
water_depth=1.0;

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

//number of cells in water
N_water=water_depth/lz[IMG+1];
printf("N_water=%d\n",N_water);

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

/*...size of the imaginary cells...*/
lx[IMG]=lx[IMG+1];
lx[IMG-1]=lx[IMG+2];
lx[IMG-2]=lx[IMG+3];
ly[IMG]=ly[IMG+1];
ly[IMG-1]=ly[IMG+2];
lz[IMG-2]=lz[IMG+3];
lz[IMG]=lz[IMG+1];
lz[IMG-1]=lz[IMG+2];
lz[IMG-2]=lz[IMG+3];
lx[NX+IMG+1]=lx[NX+IMG];
lx[NX+IMG+2]=lx[NX+IMG-1];
lx[NX+IMG+3]=lx[NX+IMG-2];
ly[NY+IMG+1]=ly[NY+IMG];
ly[NY+IMG+2]=ly[NY+IMG-1];
ly[NY+IMG+3]=ly[NY+IMG-2];
lz[NZ+IMG+1]=lz[NZ+IMG];
lz[NZ+IMG+2]=lz[NZ+IMG-1];
lz[NZ+IMG+3]=lz[NZ+IMG-2];

FILE *f1;											//file-writing pointer
FILE *fp;											//file-reading pointer for p,f-grid
FILE *fU;											//file-reading pointer for U-grid
FILE *fW;											//file-reading pointer for W-grid

/*...writing the FLUID cell-centers to file...*/
f1=fopen("CGRID.dat","w");
for (i=IMG+1;i<=NX+IMG;i++)
    for (j=IMG+1;j<=NY+IMG;j++)
        for (k=IMG+1;k<=NZ+IMG;k++)
            fprintf(f1,"%f\t%f\t%f\n",cx[i],cy[j],cz[k]);
fclose(f1);

/*...initialization of liquid and solid VOF fields...*/
for (i=0;i<=NX+IMG+5;i++)
{
    for (j=0;j<=NY+IMG+5;j++)
    {
        for (k=0;k<=NZ+IMG+5;k++)
	{
            liquid_VOF[i][j][k]=0.0;
	    solid_VOF[i][j][k]=0.0;
	}
    }
}

/*...initial patching of water upto SWL...*/
for (i=IMG+1;i<=NX+IMG;i++)
    for (j=IMG+1;j<=NY+IMG;j++)
        for (k=IMG+1;k<=N_water+IMG;k++)
            liquid_VOF[i][j][k]=1.0;

/*...reading of the solid phase - initialization of the cylinder...*/
f1=fopen("CYLINDER.dat","r");                                             			//MODIFY the FILE-NAME here
for (i=IMG+1;i<=NX+IMG;i++)
{
    for (j=IMG+1;j<=NY+IMG;j++)
    {
        for (k=IMG+1;k<=NZ+IMG;k++)
	{
            fscanf(f1,"%*f\t%*f\t%*f\t%lf\n",&solid_VOF[i][j][k]);
	    initsolid[i][j][k]=solid_VOF[i][j][k];
	}
    }
}
fclose(f1);

/*...it's either liquid or solid near and within the cylinder...*/
for (i=IMG+1;i<=NX+IMG;i++)
    for (j=IMG+1;j<=NY+IMG;j++)
        for (k=IMG+1;k<=N_water+IMG;k++)							//only need to look at the lower half of the domain
		if (solid_VOF[i][j][k]>0)
			liquid_VOF[i][j][k]=(1.0-solid_VOF[i][j][k]);

/*...hydrostatic pressure initialization below the SWL...*/					//densities of both solid and liquid need to be considered
for (i=IMG+1;i<=NX+IMG;i++)
    for (j=IMG+1;j<=NY+IMG;j++)
        for (k=IMG+1;k<=N_water+IMG;k++)
            p[i][j][k]=(liquid_VOF[i][j][k]*rho_liquid+solid_VOF[i][j][k]*rho_solid)*g*(water_depth-cz[k]);

/*...computation of the initial volume of liquid...*/
initVOL=0.0;
for (i=IMG+1;i<=NX+IMG;i++)
    for (j=IMG+1;j<=NY+IMG;j++)
        for (k=IMG+1;k<=NZ+IMG;k++)
            initVOL=initVOL+liquid_VOF[i][j][k]*lx[i]*ly[j]*lz[k];
printf("The initial volume is %0.15f\n\n",initVOL);

/*...time stepping and file-writing info...*/
tme=0.0,OLD_tme_FNPT=0.0,NEW_tme_FNPT=0.0;
double init_time_step=0.0025;									//initial time-step size for NS
dt_FNPT=0.0025;											//time-step size in the FNPT simulation - see "FEMINPUT.dat"
dt=init_time_step;										//time-step size for NS (varies according to Courant number)
int N_FNPT=5513;										//number of steps executed in the FNPT simulation - (total no. of lines in "Output_PHIT.dat")/(NX_FNPT*NZ_FNPT)
double C_PERMIT=0.2;										//max. permissible Courant number (DO NOT EXCEED 0.25)
int trigger=0;											//trigger for variable time-stepping
count=1;											//NS time-step counter
int ins=1;											//integer counter for file-writing
double tme_WRITE=5.0;										//starting time for data file-writing
int freq_WRITE=50;										//frequency of data file-writing (1 file every "freq_WRITE" time-steps)

/*...advection schemes...*/
//1-FOU (First Order Upwind) 									[unconditionally STABLE]
//2-SOU (Second Order Upwind) 									[unstable if not blended]!!!!!!!!!!!!!!!!
//3-QUICK (Quadratic Upstream Interpolation of Convection Kinetics) 				[unstable if not blended]!!!!!!!!!!!!!!!!
//4-CD2 (Central Difference) 									[unstable if not blended]!!!!!!!!!!!!!!!!
//5-TVD (Total Variation Diminishing) 								[unstable if not blended]!!!!!!!!!!!!!!!!
//6-FiOU (Fifth Order Upwind) 									[unstable if not blended]!!!!!!!!!!!!!!!!
int sch1=1;											/*...LIMITER SCHEME...*/
int sch2=6;											/*...HIGHER-ORDER SCHEME...*/
double schbln=0.3;										/*...FRACTION OF LIMITER SCHEME...*/

/******************************************************************************************************************************/
/*.......................................3D-RANS SOLVER LOOP - HARD-CODED - DO NOT CHANGE.....................................*/
/******************************************************************************************************************************/
fp=fopen("Interp_pf.dat","r");
fU=fopen("Interp_Ustag.dat","r");
fW=fopen("Interp_Wstag.dat","r");
do
{
//CALCULATION OF MAXIMUM COURANT NUMBER
double maxU=0.0,maxV=0.0,maxW=0.0,mincelX=0.0,mincelY=0.0,mincelZ=0.0;
double xcour=0.0,ycour=0.0,zcour=0.0;
double maxcour=0.0,maxvel=0.0,mincel=0.0;
for (i=IMG+1;i<=NX+IMG;i++)
{
for (j=IMG+1;j<=NY+IMG;j++)
{
for (k=IMG+1;k<=NZ+IMG;k++)
{
	if (fabs(u[i][j][k]*dt/lx[i])>xcour)							//U-momentum
	{
	xcour=fabs(u[i][j][k]*dt/lx[i]);
	maxU=fabs(u[i][j][k]);
	mincelX=lx[i];
	}

	if (fabs(v[i][j][k]*dt/ly[j])>ycour)							//V-momentum
	{
	ycour=fabs(v[i][j][k]*dt/ly[j]);
	maxV=fabs(v[i][j][k]);
	mincelY=ly[j];
	}

	if (fabs(w[i][j][k]*dt/lz[k])>zcour)							//W-momentum
	{
	zcour=fabs(w[i][j][k]*dt/lz[k]);
	maxW=fabs(w[i][j][k]);
	mincelZ=lz[k];
	}
}
}
}
if (xcour>ycour && xcour>zcour)
{
	maxcour=xcour;	maxvel=maxU;	mincel=mincelX;
}
else if (ycour>xcour && ycour>zcour)
{
	maxcour=ycour;	maxvel=maxV;	mincel=mincelY;
}
else
{
	maxcour=zcour;	maxvel=maxW;	mincel=mincelZ;
}

//VARIABLE TIME-STEPPING (dt=init_time_step unless Courant exceeds C_PERMIT)
if (trigger==0 && maxcour<C_PERMIT)								//keep equal to the initial time-step until the trigger is activated
	dt=init_time_step;
else
{
	if (maxcour<C_PERMIT)
	{
	dt=1.01*dt;										//don't directly revert back to "init_time_step" as "maxcour" may exceed "C_PERMIT"...rather ramp up "dt" slowly
		if (dt>init_time_step)
		dt=init_time_step;								//ensure that NS time-step size does not exceed dt_FNPT
	}
	else
	{
	dt=C_PERMIT*mincel/maxvel;
		if (dt>init_time_step)
		dt=init_time_step;								//ensure that NS time-step size does not exceed dt_FNPT
	}
	trigger=1;
}

//update physical time in the NS-simulation
tme=tme+dt;

//store the old FNPT field for time interpolation
if (tme>NEW_tme_FNPT)
{
	//storing pressure and volume fractions
	for (i=IMG-3;i<=IMG;i++)
	{
	for (k=IMG+1;k<=NZ+IMG;k++)
	{
		Pold[i][IMG][k]=Pnew[i][IMG][k];
		Fold[i][IMG][k]=Fnew[i][IMG][k];
	}
	}
	//storing pressure for comparison
		Pold_fnpt=Pnew_fnpt;
	//storing u-velocities
	for (i=IMG-3;i<=IMG+1;i++)
	for (k=IMG+1;k<=NZ+IMG;k++)
		Uold[i][IMG][k]=Unew[i][IMG][k];
	//storing w-velocities
	for (i=IMG-3;i<=IMG;i++)
	for (k=IMG+1;k<=NZ+IMG;k++)
		Wold[i][IMG][k]=Wnew[i][IMG][k];	    
	//storing for comparison
		Wold_fnpt=Wnew_fnpt;	
}

	/*******************************************************************************/
	/*...consistent implementation of boundary conditions for velocities and VOF...*/
	/*******************************************************************************/
	//west planar boundary [VELOCITY-INLET]
	//directly read interpolated FNPT solution to the west boundary
	if (tme>NEW_tme_FNPT)									//update u,w,p,f fields at the current FNPT time-level
	{
	printf("The code is here\n\n");
	for (i=IMG-3;i<=IMG;i++)
		for (k=IMG+1;k<=NZ+IMG;k++)
			fscanf(fp,"%*f\t%*f\t%lf\t%lf\n",&Pnew[i][IMG][k],&Fnew[i][IMG][k]);
		for (k=IMG+1;k<=NZ+IMG;k++)							//dummy loop to skip the [IMG+1] cell layer on the collocated mesh
			{
			if (k==N_water+IMG)							//cell-center just below the SWL
			fscanf(fp,"%*f\t%*f\t%lf\t%*f\n",&Pnew_fnpt);				//read pressure for comparison against the NS solution
			else
			fscanf(fp,"%*[^\n]\n");							//skip the line
			}
	for (i=IMG-3;i<=IMG+1;i++)
		for (k=IMG+1;k<=NZ+IMG;k++)	    
			fscanf(fU,"%*f\t%*f\t%lf\n",&Unew[i][IMG][k]);
	for (i=IMG-3;i<=IMG;i++)
		for (k=IMG+1;k<=NZ+IMG;k++)	    
			fscanf(fW,"%*f\t%*f\t%lf\n",&Wnew[i][IMG][k]);
		for (k=IMG+1;k<=NZ+IMG;k++)							//dummy loop to skip the [IMG+1] cell layer on the collocated mesh
			{
			if (k==N_water+IMG)							//cell-center just below the SWL
			fscanf(fW,"%*f\t%*f\t%lf\n",&Wnew_fnpt);				//read W-velocity for comparison against the NS solution
			else
			fscanf(fW,"%*[^\n]\n");							//skip the line
			}
	}
	//update physical time in the FNPT simulation IF NS-physical time exceeds the new time-level FNPT time
	//update the "times" "OLD_tme_FNPT" and "NEW_tme_FNPT" before the actual interpolation to retain consistency with the newly updated "OLD" and "NEW" FNPT data
	if (tme>NEW_tme_FNPT)
	{
		OLD_tme_FNPT=NEW_tme_FNPT;
		NEW_tme_FNPT=NEW_tme_FNPT+dt_FNPT;
	}
	/*************************************************************************************/
	/*...time-interpolation between old and new FNPT levels (ALWAYS TIME-INTERPOLATE!)...*/
	/*************************************************************************************/
	{
	if (tme>NEW_tme_FNPT || tme<OLD_tme_FNPT)
	printf("RANS dt=%f exceeds dt_FNPT=%f.....INTERPOLATION ERROR!!!!!!!!!!!!!!!\n",dt,dt_FNPT);
		//interpolating pressure and volume fractions (take care to assign volume fraction to "liquid_VOF" for FFD simulations)
		for (i=IMG-3;i<=IMG;i++)
		{
		for (k=IMG+1;k<=NZ+IMG;k++)
		{
			p[i][IMG][k]=(Pold[i][IMG][k]*(NEW_tme_FNPT-tme)+Pnew[i][IMG][k]*(tme-OLD_tme_FNPT))/dt_FNPT;
			//CRITICAL! Carefully assign to liquid_VOF[i][j][k] (NOT to f[i][j][k]) for FFD formulation-based simulations			
			liquid_VOF[i][IMG][k]=(Fold[i][IMG][k]*(NEW_tme_FNPT-tme)+Fnew[i][IMG][k]*(tme-OLD_tme_FNPT))/dt_FNPT;	
		}
		}
		//interpolating pressure for comparison
			p_fnpt=(Pold_fnpt*(NEW_tme_FNPT-tme)+Pnew_fnpt*(tme-OLD_tme_FNPT))/dt_FNPT;
		//interpolating u-velocities
		for (i=IMG-3;i<=IMG+1;i++)
		for (k=IMG+1;k<=NZ+IMG;k++)
			u[i][IMG][k]=(Uold[i][IMG][k]*(NEW_tme_FNPT-tme)+Unew[i][IMG][k]*(tme-OLD_tme_FNPT))/dt_FNPT;
		//interpolating w-velocities
		for (i=IMG-3;i<=IMG;i++)
		for (k=IMG+1;k<=NZ+IMG;k++)
			w[i][IMG][k]=(Wold[i][IMG][k]*(NEW_tme_FNPT-tme)+Wnew[i][IMG][k]*(tme-OLD_tme_FNPT))/dt_FNPT;	    
		//interpolating w-velocity for comparison
			w_fnpt=(Wold_fnpt*(NEW_tme_FNPT-tme)+Wnew_fnpt*(tme-OLD_tme_FNPT))/dt_FNPT;			
	}
        for (j=0;j<=NY+IMG+5;j++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            u[IMG+1][j][k]=u[IMG+1][IMG][k];							//stagger
            u[IMG][j][k]=u[IMG][IMG][k];							//stagger
	    u[IMG-1][j][k]=u[IMG-1][IMG][k];							//stagger
            v[IMG][j][k]=-v[IMG+1][j][k];							//collocated
            v[IMG-1][j][k]=2.0*v[IMG][j][k]-v[IMG+1][j][k];					//collocated
	    v[IMG-2][j][k]=2.0*v[IMG-1][j][k]-v[IMG][j][k];					//collocated
            w[IMG][j][k]=w[IMG][IMG][k];							//collocated
            w[IMG-1][j][k]=w[IMG-1][IMG][k];							//collocated
	    w[IMG-2][j][k]=w[IMG-2][IMG][k];							//collocated
        }
        }	
        //east planar boundary [FREE-SLIP]
        for (j=0;j<=NY+IMG+5;j++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            u[NX+IMG+1][j][k]=0.0;								//stagger
            u[NX+IMG+2][j][k]=2.0*u[NX+IMG+1][j][k]-u[NX+IMG][j][k];				//stagger
	    u[NX+IMG+3][j][k]=2.0*u[NX+IMG+2][j][k]-u[NX+IMG+1][j][k];				//stagger
            v[NX+IMG+1][j][k]=v[NX+IMG][j][k];							//collocated
            v[NX+IMG+2][j][k]=2.0*v[NX+IMG+1][j][k]-v[NX+IMG][j][k];				//collocated
	    v[NX+IMG+3][j][k]=2.0*v[NX+IMG+2][j][k]-v[NX+IMG+1][j][k];				//collocated
            w[NX+IMG+1][j][k]=w[NX+IMG][j][k];							//collocated
            w[NX+IMG+2][j][k]=2.0*w[NX+IMG+1][j][k]-w[NX+IMG][j][k];				//collocated
	    w[NX+IMG+3][j][k]=2.0*w[NX+IMG+2][j][k]-w[NX+IMG+1][j][k];				//collocated
        }
        }
        //south planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            u[i][IMG][k]=u[i][IMG+1][k];							//collocated
            u[i][IMG-1][k]=2.0*u[i][IMG][k]-u[i][IMG+1][k];					//collocated
	    u[i][IMG-2][k]=2.0*u[i][IMG-1][k]-u[i][IMG][k];					//collocated
            v[i][IMG+1][k]=0.0;									//stagger
            v[i][IMG][k]=2.0*v[i][IMG+1][k]-v[i][IMG+2][k];					//stagger
	    v[i][IMG-1][k]=2.0*v[i][IMG][k]-v[i][IMG+1][k];					//stagger
            w[i][IMG][k]=w[i][IMG+1][k];							//collocated
            w[i][IMG-1][k]=2.0*w[i][IMG][k]-w[i][IMG+1][k];					//collocated
	    w[i][IMG-2][k]=2.0*w[i][IMG-1][k]-w[i][IMG][k];					//collocated
        }
        }
        //north planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            u[i][NY+IMG+1][k]=u[i][NY+IMG][k];							//collocated
            u[i][NY+IMG+2][k]=2.0*u[i][NY+IMG+1][k]-u[i][NY+IMG][k];				//collocated
	    u[i][NY+IMG+3][k]=2.0*u[i][NY+IMG+2][k]-u[i][NY+IMG+1][k];				//collocated
            v[i][NY+IMG+1][k]=0.0;								//stagger
            v[i][NY+IMG+2][k]=2.0*v[i][NY+IMG+1][k]-v[i][NY+IMG][k];				//stagger
	    v[i][NY+IMG+3][k]=2.0*v[i][NY+IMG+2][k]-v[i][NY+IMG+1][k];				//stagger
            w[i][NY+IMG+1][k]=w[i][NY+IMG][k];							//collocated
            w[i][NY+IMG+2][k]=2.0*w[i][NY+IMG+1][k]-w[i][NY+IMG][k];				//collocated
	    w[i][NY+IMG+3][k]=2.0*w[i][NY+IMG+2][k]-w[i][NY+IMG+1][k];				//collocated
        }
        }
        //bottom planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (j=0;j<=NY+IMG+5;j++)
        {
            u[i][j][IMG]=u[i][j][IMG+1];							//collocated
            u[i][j][IMG-1]=2.0*u[i][j][IMG]-u[i][j][IMG+1];					//collocated
	    u[i][j][IMG-2]=2.0*u[i][j][IMG-1]-u[i][j][IMG];					//collocated
            v[i][j][IMG]=v[i][j][IMG+1];							//collocated
            v[i][j][IMG-1]=2.0*v[i][j][IMG]-v[i][j][IMG+1];					//collocated
	    v[i][j][IMG-2]=2.0*v[i][j][IMG-1]-v[i][j][IMG];					//collocated
            w[i][j][IMG+1]=0.0;									//stagger
            w[i][j][IMG]=2.0*w[i][j][IMG+1]-w[i][j][IMG+2];					//stagger
	    w[i][j][IMG-1]=2.0*w[i][j][IMG]-w[i][j][IMG+1];					//stagger
        }
        }
        //top planar boundary [PRESSURE OUTLET - AIR!]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (j=0;j<=NY+IMG+5;j++)
        {
            u[i][j][NZ+IMG+1]=u[i][j][NZ+IMG];							//collocated
            u[i][j][NZ+IMG+2]=2.0*u[i][j][NZ+IMG+1]-u[i][j][NZ+IMG];				//collocated
	    u[i][j][NZ+IMG+3]=2.0*u[i][j][NZ+IMG+2]-u[i][j][NZ+IMG+1];				//collocated
            v[i][j][NZ+IMG+1]=v[i][j][NZ+IMG];							//collocated
            v[i][j][NZ+IMG+2]=2.0*v[i][j][NZ+IMG+1]-v[i][j][NZ+IMG];				//collocated
	    v[i][j][NZ+IMG+3]=2.0*v[i][j][NZ+IMG+2]-v[i][j][NZ+IMG+1];				//collocated
            w[i][j][NZ+IMG+1]=w[i][j][NZ+IMG];							//stagger
            w[i][j][NZ+IMG+2]=2.0*w[i][j][NZ+IMG+1]-w[i][j][NZ+IMG];				//stagger
	    w[i][j][NZ+IMG+3]=2.0*w[i][j][NZ+IMG+2]-w[i][j][NZ+IMG+1];				//stagger 
        }
        }

	/*...cell-centered fluid properties allocation (at the Nth level)...*/
        for (k=IMG;k<=NZ+IMG+1;k++)								//domain extended by one imaginary cell at each end
        {
	for (j=IMG;j<=NY+IMG+1;j++)								//domain extended by one imaginary cell at each end
	{
	for (i=IMG;i<=NX+IMG+1;i++)								//domain extended by one imaginary cell at each end
	{
		must[i][j][k]=liquid_VOF[i][j][k]*mu_liquid+solid_VOF[i][j][k]*mu_solid+(1.0-liquid_VOF[i][j][k]-solid_VOF[i][j][k])*mu_gas;
		rhost[i][j][k]=liquid_VOF[i][j][k]*rho_liquid+solid_VOF[i][j][k]*rho_solid+(1.0-liquid_VOF[i][j][k]-solid_VOF[i][j][k])*rho_gas;
	}
	}
        }

        /********************************************************************************************/
	/*................OPERATOR-SPLIT CICSAM - first run for the liquid VOF field................*/
	/********************************************************************************************/
	CICSAM_RUN=1;
	/*...set the general VOF field "f" as being equal to the liquid VOF field...*/
	for (i=IMG-4;i<=NX+IMG;i++)
	    for (j=IMG-4;j<=NY+IMG;j++)
		for (k=IMG-4;k<=NZ+IMG;k++)
			f[i][j][k]=liquid_VOF[i][j][k];	

	/*...consistent implementation of VOF boundary conditions...*/
        //west planar boundary ["FNPT-INLET" --- no boundary condition/IMG cells required]
	//directly read interpolated FNPT solution for IMG+1, IMG+2......, IMG+5 real cells
        for (j=0;j<=NY+IMG+5;j++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            f[IMG][j][k]=f[IMG][IMG][k];							//collocated
            f[IMG-1][j][k]=f[IMG-1][IMG][k];							//collocated
        }
        }	
        //east planar boundary [FREE-SLIP]
        for (j=0;j<=NY+IMG+5;j++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            f[NX+IMG+1][j][k]=f[NX+IMG][j][k];							//collocated
            f[NX+IMG+2][j][k]=2.0*f[NX+IMG+1][j][k]-f[NX+IMG][j][k];				//collocated
        }
        }
        //south planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            f[i][IMG][k]=f[i][IMG+1][k];							//collocated
            f[i][IMG-1][k]=2.0*f[i][IMG][k]-f[i][IMG+1][k];					//collocated
        }
        }
        //north planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            f[i][NY+IMG+1][k]=f[i][NY+IMG][k];							//collocated
            f[i][NY+IMG+2][k]=2.0*f[i][NY+IMG+1][k]-f[i][NY+IMG][k];				//collocated
        }
        }
        //bottom planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (j=0;j<=NY+IMG+5;j++)
        {
            f[i][j][IMG]=f[i][j][IMG+1];							//collocated
            f[i][j][IMG-1]=2.0*f[i][j][IMG]-f[i][j][IMG+1];					//collocated
        }
        }
        //top planar boundary [PRESSURE OUTLET - AIR!]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (j=0;j<=NY+IMG+5;j++)
        {
            f[i][j][NZ+IMG+1]=-f[i][j][NZ+IMG];							//collocated
            f[i][j][NZ+IMG+2]=2.0*f[i][j][NZ+IMG+1]-f[i][j][NZ+IMG];				//collocated
        }
        }

        /*...COLOR function update...*/
        for (i=IMG+1;i<=NX+IMG;i++)
        {
            for (j=IMG+1;j<=NY+IMG;j++)
            {
                for (k=IMG+1;k<=NZ+IMG;k++)
                {
                    if (f[i][j][k]>0.5)
                        c[i][j][k]=1.0;
                    else
                        c[i][j][k]=0.0;
                }
            }
        }

        /*...ALTERNATE SWEEPING IN 3D...*/
        if ((count-1)%6==0)									//count=1,7,13,19 etc.../*...FIRST PERMUTATION...*/
        {
            //x-sweep\\..
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //y-sweep\\..
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //z-sweep\\...
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else if ((count-2)%6==0)								//count=2,8,14,20 etc.../*...SECOND PERMUTATION...*/
        {
            //z-sweep\\..
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //x-sweep\\..
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //y-sweep\\...
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else if ((count-3)%6==0)								//count=3,9,15,21 etc.../*...THIRD PERMUTATION...*/
        {
            //y-sweep\\..
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //z-sweep\\..
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //x-sweep\\...
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else if ((count-4)%6==0)								//count=4,10,16,22 etc.../*...FOURTH PERMUTATION...*/
        {
            //x-sweep\\..
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //z-sweep\\..
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //y-sweep\\...
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else if ((count-5)%6==0)								//count=5,11,17,23 etc.../*...FIFTH PERMUTATION...*/
        {
            //y-sweep\\..
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //x-sweep\\..
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //z-sweep\\...
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else 											//count=6,12,18,24 etc.../*...SIXTH PERMUTATION...*/
        {
            //z-sweep\\..
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //y-sweep\\..
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //x-sweep\\...
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }

	/*...set the liquid VOF field as being equal to the UPDATED general VOF field "f"...*/
	for (i=IMG-4;i<=NX+IMG;i++)
	    for (j=IMG-4;j<=NY+IMG;j++)
		for (k=IMG-4;k<=NZ+IMG;k++)
			liquid_VOF[i][j][k]=f[i][j][k];	

        /********************************************************************************************/
	/*................OPERATOR-SPLIT CICSAM - second run for the solid VOF field................*/
	/********************************************************************************************/
	CICSAM_RUN=2;
	/*...set the general VOF field "f" as being equal to the solid VOF field...*/
	for (i=IMG-4;i<=NX+IMG;i++)
	    for (j=IMG-4;j<=NY+IMG;j++)
		for (k=IMG-4;k<=NZ+IMG;k++)
			f[i][j][k]=solid_VOF[i][j][k];	

	/*...consistent implementation of VOF boundary conditions...*/
        //west planar boundary ["FNPT-INLET" --- no boundary condition/IMG cells required]
	//directly read interpolated FNPT solution for IMG+1, IMG+2......, IMG+5 real cells
        for (j=0;j<=NY+IMG+5;j++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            f[IMG][j][k]=f[IMG][IMG][k];							//collocated
            f[IMG-1][j][k]=f[IMG-1][IMG][k];							//collocated
        }
        }	
        //east planar boundary [FREE-SLIP]
        for (j=0;j<=NY+IMG+5;j++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            f[NX+IMG+1][j][k]=f[NX+IMG][j][k];							//collocated
            f[NX+IMG+2][j][k]=2.0*f[NX+IMG+1][j][k]-f[NX+IMG][j][k];				//collocated
        }
        }
        //south planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            f[i][IMG][k]=f[i][IMG+1][k];							//collocated
            f[i][IMG-1][k]=2.0*f[i][IMG][k]-f[i][IMG+1][k];					//collocated
        }
        }
        //north planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            f[i][NY+IMG+1][k]=f[i][NY+IMG][k];							//collocated
            f[i][NY+IMG+2][k]=2.0*f[i][NY+IMG+1][k]-f[i][NY+IMG][k];				//collocated
        }
        }
        //bottom planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (j=0;j<=NY+IMG+5;j++)
        {
            f[i][j][IMG]=f[i][j][IMG+1];							//collocated
            f[i][j][IMG-1]=2.0*f[i][j][IMG]-f[i][j][IMG+1];					//collocated
        }
        }
        //top planar boundary [PRESSURE OUTLET - AIR!]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (j=0;j<=NY+IMG+5;j++)
        {
            f[i][j][NZ+IMG+1]=-f[i][j][NZ+IMG];							//collocated
            f[i][j][NZ+IMG+2]=2.0*f[i][j][NZ+IMG+1]-f[i][j][NZ+IMG];				//collocated
        }
        }

        /*...COLOR function update...*/
        for (i=IMG+1;i<=NX+IMG;i++)
        {
            for (j=IMG+1;j<=NY+IMG;j++)
            {
                for (k=IMG+1;k<=NZ+IMG;k++)
                {
                    if (f[i][j][k]>0.5)
                        c[i][j][k]=1.0;
                    else
                        c[i][j][k]=0.0;
                }
            }
        }

        /*...ALTERNATE SWEEPING IN 3D...*/
        if ((count-1)%6==0)									//count=1,7,13,19 etc.../*...FIRST PERMUTATION...*/
        {
            //x-sweep\\..
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //y-sweep\\..
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //z-sweep\\...
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else if ((count-2)%6==0)								//count=2,8,14,20 etc.../*...SECOND PERMUTATION...*/
        {
            //z-sweep\\..
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //x-sweep\\..
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //y-sweep\\...
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else if ((count-3)%6==0)								//count=3,9,15,21 etc.../*...THIRD PERMUTATION...*/
        {
            //y-sweep\\..
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //z-sweep\\..
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //x-sweep\\...
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else if ((count-4)%6==0)								//count=4,10,16,22 etc.../*...FOURTH PERMUTATION...*/
        {
            //x-sweep\\..
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //z-sweep\\..
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //y-sweep\\...
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else if ((count-5)%6==0)								//count=5,11,17,23 etc.../*...FIFTH PERMUTATION...*/
        {
            //y-sweep\\..
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //x-sweep\\..
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //z-sweep\\...
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }
        else 											//count=6,12,18,24 etc.../*...SIXTH PERMUTATION...*/
        {
            //z-sweep\\..
            zsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //y-sweep\\..
            ysweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\..
                    retribution();

            //x-sweep\\...
            xsweep();

                //VOF and COLOR FUNCTION update\\...
                vof_c_update();

                    //redistribution\\...
                    retribution();
        }

	/*...set the solid VOF field as being equal to the UPDATED general VOF field "f"...*/
	for (i=IMG-4;i<=NX+IMG;i++)
	    for (j=IMG-4;j<=NY+IMG;j++)
		for (k=IMG-4;k<=NZ+IMG;k++)
			solid_VOF[i][j][k]=f[i][j][k];	

	/************************************************************************/
	/*.........................SEMI-EXPLICIT SOLVER.........................*/
	/************************************************************************/

	/*...Computation of X-momentum...*/
	for (i=IMG+2;i<=NX+IMG;i++)								//backward stagger
	{
	for (j=IMG+1;j<=NY+IMG;j++)
	{
        for (k=IMG+1;k<=NZ+IMG;k++)
        {
		//estimating the ADVECTING variable using a generalised CD scheme
		ue=(u[i][j][k]*(0.5*(lx[i+1]+lx[i]))+u[i+1][j][k]*(0.5*(lx[i]+lx[i-1])))/((0.5*(lx[i+1]+lx[i-1]))+lx[i]);
		uw=(u[i][j][k]*(0.5*(lx[i-1]+lx[i-2]))+u[i-1][j][k]*(0.5*(lx[i]+lx[i-1])))/((0.5*(lx[i]+lx[i-2]))+lx[i-1]);
		un=(u[i][j][k]*ly[j+1]+u[i][j+1][k]*ly[j])/(ly[j]+ly[j+1]);
		us=(u[i][j][k]*ly[j-1]+u[i][j-1][k]*ly[j])/(ly[j]+ly[j-1]);
		ut=(u[i][j][k]*lz[k+1]+u[i][j][k+1]*lz[k])/(lz[k]+lz[k+1]);
		ub=(u[i][j][k]*lz[k-1]+u[i][j][k-1]*lz[k])/(lz[k]+lz[k-1]);
		vn=(v[i][j+1][k]*lx[i-1]+v[i-1][j+1][k]*lx[i])/(lx[i-1]+lx[i]);			//stagger
		vs=(v[i][j][k]*lx[i-1]+v[i-1][j][k]*lx[i])/(lx[i-1]+lx[i]);			//stagger
		wt=(w[i][j][k+1]*lx[i-1]+w[i-1][j][k+1]*lx[i])/(lx[i-1]+lx[i]);			//stagger
		wb=(w[i][j][k]*lx[i-1]+w[i-1][j][k]*lx[i])/(lx[i-1]+lx[i]);			//stagger

		//advection of x-momentum
		sch=sch1;									//estimating advected momentum using the LIMITER scheme
		uadrou();
		FXAD1=(dt/(0.5*(lx[i]+lx[i-1])))*(uwa*uw-uea*ue)+(dt/ly[j])*(usa*vs-una*vn)+(dt/lz[k])*(uba*wb-uta*wt);
		sch=sch2;									//estimating advected momentum using the HIGHER-ORDER scheme
		uadrou();
		FXAD2=(dt/(0.5*(lx[i]+lx[i-1])))*(uwa*uw-uea*ue)+(dt/ly[j])*(usa*vs-una*vn)+(dt/lz[k])*(uba*wb-uta*wt);
		//blending of LIMITER and HIGHER-ORDER schemes
		FXAD=schbln*FXAD1+(1.0-schbln)*FXAD2;

		//evaluation of density at U-momentum cell centre
		rhow=(rhost[i][j][k]*lx[i-1]+rhost[i-1][j][k]*lx[i])/(lx[i-1]+lx[i]);

		//two-phase viscosity [weighted Harmonic averaging]
	        mustar_FX();

		//evaluation of the U-momentum cell center distances
		dxe=lx[i];									//stagger
		dxw=lx[i-1];									//stagger
		dyn=0.5*(ly[j]+ly[j+1]);
		dys=0.5*(ly[j]+ly[j-1]);
		dzt=0.5*(lz[k]+lz[k+1]);
		dzb=0.5*(lz[k]+lz[k-1]);

		//diffusion of x-momentum
		FXDF=((dt*mue)/rhow)*((u[i+1][j][k]-u[i][j][k])/(dxe*(0.5*(lx[i]+lx[i-1]))))+((dt*muw)/rhow)*((u[i-1][j][k]-u[i][j][k])/(dxw*(0.5*(lx[i]+lx[i-1]))))+((dt*mun)/rhow)*((u[i][j+1][k]-u[i][j][k])/(dyn*ly[j]))+((dt*mus)/rhow)*((u[i][j-1][k]-u[i][j][k])/(dys*ly[j]))+((dt*mut)/rhow)*((u[i][j][k+1]-u[i][j][k])/(dzt*lz[k]))+((dt*mub)/rhow)*((u[i][j][k-1]-u[i][j][k])/(dzb*lz[k]));
		//EVALUATION OF MOMENTUM CELL CENTRED TILDE VELOCITIES
		utild[i][j][k]=u[i][j][k]+FXAD+FXDF;
	}
	}
	}

	/*...Computation of Y-momentum...*/
	for (i=IMG+1;i<=NX+IMG;i++)
	{
	for (j=IMG+2;j<=NY+IMG;j++)								//backward stagger
	{
        for (k=IMG+1;k<=NZ+IMG;k++)
        {
		//estimating the ADVECTING variable using a generalised CD scheme
		ue=(u[i+1][j][k]*ly[j-1]+u[i+1][j-1][k]*ly[j])/(ly[j-1]+ly[j]);			//stagger
		uw=(u[i][j][k]*ly[j-1]+u[i][j-1][k]*ly[j])/(ly[j-1]+ly[j]);			//stagger
		ve=(v[i][j][k]*lx[i+1]+v[i+1][j][k]*lx[i])/(lx[i+1]+lx[i]);
		vw=(v[i][j][k]*lx[i-1]+v[i-1][j][k]*lx[i])/(lx[i-1]+lx[i]);
		vn=(v[i][j][k]*(0.5*(ly[j+1]+ly[j]))+v[i][j+1][k]*(0.5*(ly[j]+ly[j-1])))/((0.5*(ly[j+1]+ly[j-1]))+ly[j]);
		vs=(v[i][j][k]*(0.5*(ly[j-1]+ly[j-2]))+v[i][j-1][k]*(0.5*(ly[j]+ly[j-1])))/((0.5*(ly[j]+ly[j-2]))+ly[j-1]);
		vt=(v[i][j][k]*lz[k+1]+v[i][j][k+1]*lz[k])/(lz[k+1]+lz[k]);
		vb=(v[i][j][k]*lz[k-1]+v[i][j][k-1]*lz[k])/(lz[k-1]+lz[k]);
                wt=(w[i][j][k+1]*ly[j-1]+w[i][j-1][k+1]*ly[j])/(ly[j-1]+ly[j]);			//stagger
                wb=(w[i][j][k]*ly[j-1]+w[i][j-1][k]*ly[j])/(ly[j-1]+ly[j]);			//stagger

		//advection of y-momentum
		sch=sch1;									//estimating advected momentum using the LIMITER scheme
		vadrou();
		FYAD1=(dt/lx[i])*(vwa*uw-vea*ue)+(dt/(0.5*(ly[j]+ly[j-1])))*(vsa*vs-vna*vn)+(dt/lz[k])*(vba*wb-vta*wt);
		sch=sch2;									//estimating advected momentum using the HIGHER-ORDER scheme
		vadrou();
		FYAD2=(dt/lx[i])*(vwa*uw-vea*ue)+(dt/(0.5*(ly[j]+ly[j-1])))*(vsa*vs-vna*vn)+(dt/lz[k])*(vba*wb-vta*wt);
		//blending of LIMITER and HIGHER-ORDER schemes
		FYAD=schbln*FYAD1+(1.0-schbln)*FYAD2;

		//evaluation of density at V-momentum cell centre
		rhos=(rhost[i][j][k]*ly[j-1]+rhost[i][j-1][k]*ly[j])/(ly[j-1]+ly[j]);

		//two-phase viscosity [weighted Harmonic averaging]
                mustar_FY();

		//evaluation of the momentum cell centre distances
		dxe=0.5*(lx[i]+lx[i+1]);
		dxw=0.5*(lx[i]+lx[i-1]);
		dyn=ly[j];									//stagger
		dys=ly[j-1];									//stagger
                dzt=0.5*(lz[k]+lz[k+1]);
                dzb=0.5*(lz[k]+lz[k-1]);

		//diffusion of y-momentum
		FYDF=((dt*mue)/rhos)*((v[i+1][j][k]-v[i][j][k])/(dxe*lx[i]))+((dt*muw)/rhos)*((v[i-1][j][k]-v[i][j][k])/(dxw*lx[i]))+((dt*mun)/rhos)*((v[i][j+1][k]-v[i][j][k])/(dyn*(0.5*(ly[j]+ly[j-1]))))+((dt*mus)/rhos)*((v[i][j-1][k]-v[i][j][k])/(dys*(0.5*(ly[j]+ly[j-1]))))+((dt*mut)/rhos)*((v[i][j][k+1]-v[i][j][k])/(dzt*lz[k]))+((dt*mub)/rhos)*((v[i][j][k-1]-v[i][j][k])/(dzb*lz[k]));
		//EVALUATION OF MOMENTUM CELL CENTRED TILDE VELOCITIES
		vtild[i][j][k]=v[i][j][k]+FYAD+FYDF;
	}
	}
	}

	/*...Computation of Z-momentum...*/
        for (i=IMG+1;i<=NX+IMG;i++)
        {
        for (j=IMG+1;j<=NY+IMG;j++)
        {
        for (k=IMG+2;k<=NZ+IMG;k++)								//backward stagger
        {
		//estimating the ADVECTING variable using a generalised CD scheme
		ue=(u[i+1][j][k]*lz[k-1]+u[i+1][j][k-1]*lz[k])/(lz[k-1]+lz[k]);			//stagger
		uw=(u[i][j][k]*lz[k-1]+u[i][j][k-1]*lz[k])/(lz[k-1]+lz[k]);			//stagger
		vn=(v[i][j+1][k]*lz[k-1]+v[i][j+1][k-1]*lz[k])/(lz[k-1]+lz[k]);			//stagger
		vs=(v[i][j][k]*lz[k-1]+v[i][j][k-1]*lz[k])/(lz[k-1]+lz[k]);			//stagger
		we=(w[i][j][k]*lx[i+1]+w[i+1][j][k]*lx[i])/(lx[i+1]+lx[i]);
		ww=(w[i][j][k]*lx[i-1]+w[i-1][j][k]*lx[i])/(lx[i-1]+lx[i]);
		wn=(w[i][j][k]*ly[j+1]+w[i][j+1][k]*ly[j])/(ly[j+1]+ly[j]);
		ws=(w[i][j][k]*ly[j-1]+w[i][j-1][k]*ly[j])/(ly[j-1]+ly[j]);
		wt=(w[i][j][k]*(0.5*(lz[k+1]+lz[k]))+w[i][j][k+1]*(0.5*(lz[k]+lz[k-1])))/((0.5*(lz[k+1]+lz[k-1]))+lz[k]);
		wb=(w[i][j][k]*(0.5*(lz[k-1]+lz[k-2]))+w[i][j][k-1]*(0.5*(lz[k]+lz[k-1])))/((0.5*(lz[k]+lz[k-2]))+lz[k-1]);

        	//advection of z-momentum
		sch=sch1;									//estimating advected momentum using the LIMITER scheme
		wadrou();
		FZAD1=(dt/lx[i])*(wwa*uw-wea*ue)+(dt/ly[j])*(wsa*vs-wna*vn)+(dt/(0.5*(lz[k]+lz[k-1])))*(wba*wb-wta*wt);
		sch=sch2;									//estimating advected momentum using the HIGHER-ORDER scheme
		wadrou();
		FZAD2=(dt/lx[i])*(wwa*uw-wea*ue)+(dt/ly[j])*(wsa*vs-wna*vn)+(dt/(0.5*(lz[k]+lz[k-1])))*(wba*wb-wta*wt);
		//blending of LIMITER and HIGHER-ORDER schemes
		FZAD=schbln*FZAD1+(1.0-schbln)*FZAD2;

		//evaluation of density at W-momentum cell centre
		rhob=(rhost[i][j][k]*lz[k-1]+rhost[i][j][k-1]*lz[k])/(lz[k-1]+lz[k]);

		//two-phase viscosity [weighted Harmonic averaging]
                mustar_FZ();

		//evaluation of the momentum cell centre distances
		dxe=0.5*(lx[i]+lx[i+1]);
		dxw=0.5*(lx[i]+lx[i-1]);
		dyn=0.5*(ly[j]+ly[j+1]);
		dys=0.5*(ly[j]+ly[j-1]);
                dzt=lz[k];									//stagger
                dzb=lz[k-1];									//stagger

		//diffusion of z-momentum
		FZDF=((dt*mue)/rhob)*((w[i+1][j][k]-w[i][j][k])/(dxe*lx[i]))+((dt*muw)/rhob)*((w[i-1][j][k]-w[i][j][k])/(dxw*lx[i]))+((dt*mun)/rhob)*((w[i][j+1][k]-w[i][j][k])/(dyn*ly[j]))+((dt*mus)/rhob)*((w[i][j-1][k]-w[i][j][k])/(dys*ly[j]))+((dt*mut)/rhob)*((w[i][j][k+1]-w[i][j][k])/(dzt*(0.5*(lz[k]+lz[k-1]))))+((dt*mub)/rhob)*((w[i][j][k-1]-w[i][j][k])/(dzb*(0.5*(lz[k]+lz[k-1]))));
		//source of z-momentum due to gravity
		FZSRC_GR=-g*dt;
		//EVALUATION OF MOMENTUM CELL CENTRED TILDE VELOCITIES
		wtild[i][j][k]=w[i][j][k]+FZAD+FZDF+FZSRC_GR;
        }
        }
        }

	/*...consistent implementation of pressure boundary conditions...*/
        //west planar boundary [INLET]
        for (j=0;j<=NY+IMG+5;j++)
            for (k=0;k<=NZ+IMG+5;k++)
                p[IMG][j][k]=p[IMG][IMG][k];							//collocated
        //east planar boundary [FREE-SLIP]
        for (j=0;j<=NY+IMG+5;j++)
            for (k=0;k<=NZ+IMG+5;k++)
                p[NX+IMG+1][j][k]=p[NX+IMG][j][k];						//zero-gradient
        //south planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
            for (k=0;k<=NZ+IMG+5;k++)
                p[i][IMG][k]=p[i][IMG+1][k];							//zero-gradient
        //north planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
            for (k=0;k<=NZ+IMG+5;k++)
                p[i][NY+IMG+1][k]=p[i][NY+IMG][k];						//zero-gradient
        //bottom planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
            for (j=0;j<=NY+IMG+5;j++)
                p[i][j][IMG]=p[i][j][IMG+1];							//zero-gradient
        //top planar boundary [PRESSURE OUTLET - AIR!]
        for (i=0;i<=NX+IMG+5;i++)
            for (j=0;j<=NY+IMG+5;j++)
                p[i][j][NZ+IMG+1]=-p[i][j][NZ+IMG];						//zero-gauge pressure

	/*...obtaining x-direction starred (*) velocities...*/
        for (k=IMG+1;k<=NZ+IMG;k++)
        {
	for (j=IMG+1;j<=NY+IMG;j++)
	{
	for (i=IMG+2;i<=NX+IMG;i++)								//backward stagger
	{
            rhow=(rhost[i][j][k]*lx[i-1]+rhost[i-1][j][k]*lx[i])/(lx[i-1]+lx[i]);
            ustar[i][j][k]=utild[i][j][k]-(dt/(0.5*rhow*(lx[i]+lx[i-1])))*(p[i][j][k]-p[i-1][j][k]);//stagger
	}
	}
        }

	/*...obtaining y-direction starred (*) velocities...*/
        for (k=IMG+1;k<=NZ+IMG;k++)
        {
	for (j=IMG+2;j<=NY+IMG;j++)								//backward stagger
	{
	for (i=IMG+1;i<=NX+IMG;i++)
	{
            rhos=(rhost[i][j][k]*ly[j-1]+rhost[i][j-1][k]*ly[j])/(ly[j-1]+ly[j]);
            vstar[i][j][k]=vtild[i][j][k]-(dt/(0.5*rhos*(ly[j]+ly[j-1])))*(p[i][j][k]-p[i][j-1][k]);//stagger
	}
	}
        }

        /*...obtaining z-direction starred (*) velocities...*/
        for (k=IMG+2;k<=NZ+IMG;k++)								//backward stagger
        {
	for (j=IMG+1;j<=NY+IMG;j++)
	{
	for (i=IMG+1;i<=NX+IMG;i++)
	{
            rhob=(rhost[i][j][k]*lz[k-1]+rhost[i][j][k-1]*lz[k])/(lz[k-1]+lz[k]);
            wstar[i][j][k]=wtild[i][j][k]-(dt/(0.5*rhob*(lz[k]+lz[k-1])))*(p[i][j][k]-p[i][j][k-1]);//stagger
	}
	}
        }

	/*...consistent implementation of starred (*) velocity boundary conditions...*/
	//west planar boundary [VELOCITY-INLET]
        for (j=0;j<=NY+IMG+5;j++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            ustar[IMG+1][j][k]=u[IMG+1][IMG][k];						//stagger
            ustar[IMG][j][k]=u[IMG][IMG][k];							//stagger
            vstar[IMG][j][k]=-vstar[IMG+1][j][k];						//collocated
            vstar[IMG-1][j][k]=2.0*vstar[IMG][j][k]-vstar[IMG+1][j][k];				//collocated
            wstar[IMG][j][k]=w[IMG][IMG][k];							//collocated
            wstar[IMG-1][j][k]=w[IMG-1][IMG][k];						//collocated
        }
        }	
        //east planar boundary [FREE-SLIP]
        for (j=0;j<=NY+IMG+5;j++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            ustar[NX+IMG+1][j][k]=0.0;								//stagger
            ustar[NX+IMG+2][j][k]=2.0*ustar[NX+IMG+1][j][k]-ustar[NX+IMG][j][k];		//stagger
            vstar[NX+IMG+1][j][k]=vstar[NX+IMG][j][k];						//collocated
            vstar[NX+IMG+2][j][k]=2.0*vstar[NX+IMG+1][j][k]-vstar[NX+IMG][j][k];		//collocated
            wstar[NX+IMG+1][j][k]=wstar[NX+IMG][j][k];						//collocated
            wstar[NX+IMG+2][j][k]=2.0*wstar[NX+IMG+1][j][k]-wstar[NX+IMG][j][k];		//collocated
        }
        }
        //south planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            ustar[i][IMG][k]=ustar[i][IMG+1][k];						//collocated
            ustar[i][IMG-1][k]=2.0*ustar[i][IMG][k]-ustar[i][IMG+1][k];				//collocated
            vstar[i][IMG+1][k]=0.0;								//stagger
            vstar[i][IMG][k]=2.0*vstar[i][IMG+1][k]-vstar[i][IMG+2][k];				//stagger
            wstar[i][IMG][k]=wstar[i][IMG+1][k];						//collocated
            wstar[i][IMG-1][k]=2.0*wstar[i][IMG][k]-wstar[i][IMG+1][k];				//collocated
        }
        }
        //north planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (k=0;k<=NZ+IMG+5;k++)
        {
            ustar[i][NY+IMG+1][k]=ustar[i][NY+IMG][k];						//collocated
            ustar[i][NY+IMG+2][k]=2.0*ustar[i][NY+IMG+1][k]-ustar[i][NY+IMG][k];		//collocated
            vstar[i][NY+IMG+1][k]=0.0;								//stagger
            vstar[i][NY+IMG+2][k]=2.0*vstar[i][NY+IMG+1][k]-vstar[i][NY+IMG][k];		//stagger
            wstar[i][NY+IMG+1][k]=wstar[i][NY+IMG][k];						//collocated
            wstar[i][NY+IMG+2][k]=2.0*wstar[i][NY+IMG+1][k]-wstar[i][NY+IMG][k];		//collocated
        }
        }
        //bottom planar boundary [FREE-SLIP]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (j=0;j<=NY+IMG+5;j++)
        {
            ustar[i][j][IMG]=ustar[i][j][IMG+1];						//collocated
            ustar[i][j][IMG-1]=2.0*ustar[i][j][IMG]-ustar[i][j][IMG+1];				//collocated
            vstar[i][j][IMG]=vstar[i][j][IMG+1];						//collocated
            vstar[i][j][IMG-1]=2.0*vstar[i][j][IMG]-vstar[i][j][IMG+1];				//collocated
            wstar[i][j][IMG+1]=0.0;								//stagger
            wstar[i][j][IMG]=2.0*wstar[i][j][IMG+1]-wstar[i][j][IMG+2];				//stagger
        }
        }
        //top planar boundary [PRESSURE OUTLET - AIR!]
        for (i=0;i<=NX+IMG+5;i++)
        {
        for (j=0;j<=NY+IMG+5;j++)
        {
            ustar[i][j][NZ+IMG+1]=ustar[i][j][NZ+IMG];						//collocated
            ustar[i][j][NZ+IMG+2]=2.0*ustar[i][j][NZ+IMG+1]-ustar[i][j][NZ+IMG];		//collocated
            vstar[i][j][NZ+IMG+1]=vstar[i][j][NZ+IMG];						//collocated
            vstar[i][j][NZ+IMG+2]=2.0*vstar[i][j][NZ+IMG+1]-vstar[i][j][NZ+IMG];		//collocated
            wstar[i][j][NZ+IMG+1]=wstar[i][j][NZ+IMG];						//stagger
            wstar[i][j][NZ+IMG+2]=2.0*wstar[i][j][NZ+IMG+1]-wstar[i][j][NZ+IMG];		//stagger
        }
        }

	/*...obtaining the MASS SOURCE for all cells...*/
        for (k=IMG+1;k<=NZ+IMG;k++)
	for (j=IMG+1;j<=NY+IMG;j++)
	for (i=IMG+1;i<=NX+IMG;i++)
		mass[i][j][k]=(ustar[i][j][k]-ustar[i+1][j][k])*ly[j]*lz[k]+(vstar[i][j][k]-vstar[i][j+1][k])*lx[i]*lz[k]+(wstar[i][j][k]-wstar[i][j][k+1])*lx[i]*ly[j];//backward stagger

	/*******************************************************************************/
	/*..............the equation of pressure correction (EOPC).....................*/
	/*******************************************************************************/
	/*...initializing all corrections from zero...*/
        for (k=IMG+1;k<=NZ+IMG;k++)
        for (j=IMG+1;j<=NY+IMG;j++)
        for (i=IMG+1;i<=NX+IMG;i++)
		pcorr[i][j][k]=0.0;

	gcount=1;
	/*the GSSOR solution loop*/
	printf("\n..............Equation of pressure correction begins..........................\n");
	do
	{
  	    gsum=0.0;
	    /*...consistent implementation of pressure correction boundary conditions...*/
            //west planar boundary [INLET]
            for (j=0;j<=NY+IMG+5;j++)
                for (k=0;k<=NZ+IMG+5;k++)
                    pcorr[IMG][j][k]=pcorr[IMG+1][j][k];					//zero-gradient
            //east planar boundary [FREE-SLIP]
            for (j=0;j<=NY+IMG+5;j++)
                for (k=0;k<=NZ+IMG+5;k++)
                    pcorr[NX+IMG+1][j][k]=pcorr[NX+IMG][j][k];					//zero-gradient
            //south planar boundary [FREE-SLIP]
            for (i=0;i<=NX+IMG+5;i++)
                for (k=0;k<=NZ+IMG+5;k++)
                    pcorr[i][IMG][k]=pcorr[i][IMG+1][k];					//zero-gradient
            //north planar boundary [FREE-SLIP]
            for (i=0;i<=NX+IMG+5;i++)
                for (k=0;k<=NZ+IMG+5;k++)
                    pcorr[i][NY+IMG+1][k]=pcorr[i][NY+IMG][k];					//zero-gradient
            //bottom planar boundary [FREE-SLIP]
            for (i=0;i<=NX+IMG+5;i++)
                for (j=0;j<=NY+IMG+5;j++)
                    pcorr[i][j][IMG]=pcorr[i][j][IMG+1];					//zero-gradient
            //top planar boundary [PRESSURE OUTLET - AIR!]
            for (i=0;i<=NX+IMG+5;i++)
                for (j=0;j<=NY+IMG+5;j++)
                    pcorr[i][j][NZ+IMG+1]=-pcorr[i][j][NZ+IMG];					//zero-gauge pressure

	    /*...solution of the EOPC...*/
	    for (k=IMG+1;k<=NZ+IMG;k++)
	    {
	    for (j=IMG+1;j<=NY+IMG;j++)
	    {
	    for (i=IMG+1;i<=NX+IMG;i++)
	    {
		dxe=0.5*(lx[i]+lx[i+1]);
		dxw=0.5*(lx[i]+lx[i-1]);
		dyn=0.5*(ly[j]+ly[j+1]);
		dys=0.5*(ly[j]+ly[j-1]);
		dzt=0.5*(lz[k]+lz[k+1]);
		dzb=0.5*(lz[k]+lz[k-1]);
		rhoe=(rhost[i][j][k]*lx[i+1]+rhost[i+1][j][k]*lx[i])/(lx[i+1]+lx[i]);
		rhow=(rhost[i][j][k]*lx[i-1]+rhost[i-1][j][k]*lx[i])/(lx[i-1]+lx[i]);
		rhon=(rhost[i][j][k]*ly[j+1]+rhost[i][j+1][k]*ly[j])/(ly[j+1]+ly[j]);
		rhos=(rhost[i][j][k]*ly[j-1]+rhost[i][j-1][k]*ly[j])/(ly[j-1]+ly[j]);
		rhot=(rhost[i][j][k]*lz[k+1]+rhost[i][j][k+1]*lz[k])/(lz[k+1]+lz[k]);
		rhob=(rhost[i][j][k]*lz[k-1]+rhost[i][j][k-1]*lz[k])/(lz[k-1]+lz[k]);
		apd=(dt*ly[j]*lz[k])*((1.0/(dxe*rhoe))+(1.0/(dxw*rhow)))+(dt*lx[i]*lz[k])*((1.0/(dyn*rhon))+(1.0/(dys*rhos)))+(dt*lx[i]*ly[j])*((1.0/(dzt*rhot))+(1.0/(dzb*rhob)));
		res[i][j][k]=-apd*pcorr[i][j][k]+(mass[i][j][k]+((pcorr[i+1][j][k]/(dxe*rhoe))+(pcorr[i-1][j][k]/(dxw*rhow)))*(dt*ly[j]*lz[k])+((pcorr[i][j+1][k]/(dyn*rhon))+(pcorr[i][j-1][k]/(dys*rhos)))*(dt*lx[i]*lz[k])+((pcorr[i][j][k+1]/(dzt*rhot))+(pcorr[i][j][k-1]/(dzb*rhob)))*(dt*lx[i]*ly[j]));
		corr=res[i][j][k]/(apd*(1.0+ap));
		gsq=res[i][j][k]*res[i][j][k];
		gsum=gsum+gsq;
		pcorr[i][j][k]=pcorr[i][j][k]+corr;
	    }
	    }
	    }
	    gRMS=sqrt(gsum/(NX*NY*NZ));
	    gcount=gcount+1;
	    printf ("time=%0.4f\tpcount=%d\tgRMS=%0.10f\tuRMS=%0.8f\tvRMS=%0.8f\twRMS=%0.8f\tmaxCOUR=%0.4f\tU_IMG+1=%0.4f\tU_IMG+2=%0.4f\n",tme,gcount,gRMS,uRMS,vRMS,wRMS,maxcour,u[IMG+1][IMG+1][IMG+(NZ/2)],u[IMG+2][IMG+1][IMG+(NZ/2)]);		
	}while (gRMS>EOPC_ERROR);
	printf("\n..............Equation of pressure correction ends..........................\n");

	/*...obtaining the (n+1)th level pressures...*/
	for (i=IMG+1;i<=NX+IMG;i++)
            for (j=IMG+1;j<=NY+IMG;j++)
                for (k=IMG+1;k<=NZ+IMG;k++)
                    p[i][j][k]=p[i][j][k]+pcorr[i][j][k];

	/*...obtaining the velocities at (n+1)th level...*/
	usum=0.0;
	vsum=0.0;
	wsum=0.0;

	/*...obtaining (n+1)th level U-velocities...*/
	for (i=IMG+2;i<=NX+IMG;i++)								//backward stagger
	{
	for (j=IMG+1;j<=NY+IMG;j++)
	{
        for (k=IMG+1;k<=NZ+IMG;k++)
        {
	    OLD_U=u[i][j][k];
            rhow=(rhost[i][j][k]*lx[i-1]+rhost[i-1][j][k]*lx[i])/(lx[i-1]+lx[i]);
            u[i][j][k]=utild[i][j][k]-(dt/(0.5*rhow*(lx[i]+lx[i-1])))*(p[i][j][k]-p[i-1][j][k]);//stagger
            udiff=u[i][j][k]-OLD_U;
            u[i][j][k]=u[i][j][k]*(1.0-solid_VOF[i][j][k]);					//resetting velocities within the SOLID to zero
            usq=udiff*udiff;
            usum=usum+usq;
	}
	}
	}

	/*...obtaining (n+1)th level V-velocities...*/
	for (i=IMG+1;i<=NX+IMG;i++)
	{
	for (j=IMG+2;j<=NY+IMG;j++)								//backward stagger
	{
        for (k=IMG+1;k<=NZ+IMG;k++)
        {
	    OLD_V=v[i][j][k];
            rhos=(rhost[i][j][k]*ly[j-1]+rhost[i][j-1][k]*ly[j])/(ly[j-1]+ly[j]);
            v[i][j][k]=vtild[i][j][k]-(dt/(0.5*rhos*(ly[j]+ly[j-1])))*(p[i][j][k]-p[i][j-1][k]);//stagger
            vdiff=v[i][j][k]-OLD_V;
            v[i][j][k]=v[i][j][k]*(1.0-solid_VOF[i][j][k]);					//resetting velocities within the SOLID to zero
            vsq=vdiff*vdiff;
            vsum=vsum+vsq;
	}
	}
	}

	/*...obtaining (n+1)th level W-velocities...*/
	for (i=IMG+1;i<=NX+IMG;i++)
	{
	for (j=IMG+1;j<=NY+IMG;j++)
	{
        for (k=IMG+2;k<=NZ+IMG;k++)								//backward stagger
        {
	    OLD_W=w[i][j][k];
            rhob=(rhost[i][j][k]*lz[k-1]+rhost[i][j][k-1]*lz[k])/(lz[k-1]+lz[k]);
            w[i][j][k]=wtild[i][j][k]-(dt/(0.5*rhob*(lz[k]+lz[k-1])))*(p[i][j][k]-p[i][j][k-1]);//stagger
            wdiff=w[i][j][k]-OLD_W;
            w[i][j][k]=w[i][j][k]*(1.0-solid_VOF[i][j][k]);					//resetting velocities within the SOLID to zero
            wsq=wdiff*wdiff;
            wsum=wsum+wsq;
	}
	}
	}

uRMS=sqrt(usum/((NX-1)*NY*NZ));
vRMS=sqrt(vsum/(NX*(NY-1)*NZ));
wRMS=sqrt(wsum/(NX*NY*(NZ-1)));

//computation of the final volume of fluid
finVOL=0.0;
for (i=IMG+1;i<=NX+IMG;i++)
    for (j=IMG+1;j<=NY+IMG;j++)
        for (k=IMG+1;k<=NZ+IMG;k++)
            finVOL=finVOL+liquid_VOF[i][j][k]*lx[i]*ly[j]*lz[k];

/*...Writing the VE to file...*/
f1=fopen("VE_SOLITON.dat","a");
fprintf(f1,"%f\t%0.15f\n",tme,((initVOL-finVOL)/(initVOL))*100.0);
fclose(f1);

/*...Comparison of the FNPT and NS solutions at the first real cell...*/
f1=fopen("FNPT_vs_NS_p.dat","a");
fprintf(f1,"%f\t%f\t%f\n",tme,p_fnpt,p[IMG+1][IMG+(NY/2)][IMG+N_water]);
fclose(f1);

f1=fopen("FNPT_vs_NS_w.dat","a");
fprintf(f1,"%f\t%f\t%f\n",tme,w_fnpt,w[IMG+1][IMG+(NY/2)][IMG+N_water]);
fclose(f1);

//printf ("time=%0.4f\tpcount=%d\tgRMS=%0.10f\tuRMS=%0.8f\tvRMS=%0.8f\twRMS=%0.8f\tmaxCOUR=%0.4f\tdt=%f\n",tme,gcount,gRMS,uRMS,vRMS,wRMS,maxcour,dt);		

/********************/
/*...FILE WRITING...*/
/********************/
if ((tme>tme_WRITE && count%freq_WRITE==0))									//begin FILE-WRITING
{
	sprintf (string, "000_%d.dat",ins);
	f1=fopen(string,"w");
	fprintf (f1,"\nVARIABLES = ""X"", ""Y"", ""Z"", ""ZED"" ""liquid_VOF"", ""solid_VOF""\n");
	fprintf (f1,"ZONE I=%d, J=%d, K=%d, ZONETYPE=ORDERED, DATAPACKING=POINT\n",NZ,NY,NX);
	for (i=IMG+1;i<=NX+IMG;i++)
	for (j=IMG+1;j<=NY+IMG;j++)
    	for (k=IMG+1;k<=NZ+IMG;k++)
		fprintf(f1,"%f\t%f\t%f\t%f\t%f\t%f\n",cx[i],cy[j],cz[k],cz[k],liquid_VOF[i][j][k],initsolid[i][j][k]);
	fclose(f1);
ins=ins+1;
}												//end FILE-WRITING
count=count+1;
}while(tme<dt_FNPT*N_FNPT);									//end time loop
fclose(fp);
fclose(fU);
fclose(fW);
tt=clock()-tt;
printf("Percentage mass-loss is %0.12f\n",((initVOL-finVOL)/(initVOL))*100.0);
printf("The computation time is %f seconds\n",(float)tt/CLOCKS_PER_SEC);
printf("Total number of (a)undershoot cases=%d and (b)overshoot cases=%d\n",usc,osc);
printf("Total number of redistribution cases=%d\n",usc+osc);
}												//main() ends here


