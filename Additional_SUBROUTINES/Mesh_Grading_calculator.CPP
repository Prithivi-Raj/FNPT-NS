//THIS CODE SOLVES FOR THE CELL-COUNT (nx) AND STRETCHING-FACTOR (rx) REQUIRED TO GENERATE A NON-UNIFORM MESH OVER A DOMAIN OF LENGTH "L"
//THE CODE IS ONLY INTENDED FOR A TYPICAL MESH GENERATION SCENARIO ENCOUNTERED WHERE ONLY THE DOMAIN-LENGTH, STARTING AND ENDING CELL-SIZES ARE KNOWN  
//FOLLOWING GEOMETRIC SERIES IS CONSIDERED: (1 + rx + rx^2 + rx^3 + .... + rx^(nx-1))*deltax_start = L ; where deltax_end = deltax_start*(rx^(nx-1)) 
//USE THE "-O2" OPTIMIZATION FLAG FOR FASTER EXECUTION

#include <stdio.h>
#include <math.h>

main()
{
//knowns
double L=4.0;											//length of the domain
double deltax_start=0.01;									//starting cell-size
double deltax_end=0.04;										//ending cell-size

//unknowns
double rx;											//mesh-stretching factor [UNKNOWN]
int nx;												//number of cells [UNKNOWN]

//bounds of cell-count
int nx_lower,nx_upper;										
if (deltax_start>deltax_end)									//if mesh is to be successively refined
	{
	nx_lower=L/deltax_start;								//uniformly coarse mesh throughout
	nx_upper=L/deltax_end;									//uniformly fine mesh throughout
	}
else												//if mesh is to be successively coarsened
	{
	nx_lower=L/deltax_end;									//uniformly coarse mesh throughout
	nx_upper=L/deltax_start;								//uniformly fine mesh throughout
	}
 
double tolerance=1e-03;										//error tolerance desired in the final solution
double mag_drx=1e-05;										//magnitude by which "rx" is incremented/decremented

//iterative solution routine for rx and nx
for (int count=nx_lower+1 ; count<=nx_upper-1 ; count++)					//loop between the lower and upper bounds for "nx"; excluding both values
{
printf("\nThe code is here..............count=%d............................\n",count);
rx=1.0;												// "rx" begins from unity and is incremented/decremented depending on the problem at hand
double drx=((deltax_end-deltax_start)/fabs(deltax_end-deltax_start))*(mag_drx);			//decides whether "rx" is to be incremented or decremented
double diff1=1.0,diff2=1.0;
double min1=1.0,min2=1.0,rxmin1,rxmin2;
double exsum;
	do
	{
	exsum=1.0;										//initiation of geometric progression
		diff1=deltax_start*pow(rx,count-1)-deltax_end;					//CONSTRAINT 1: relationship between first and last cell sizes
		{										//store the minimum value of "diff1"
			if (fabs(min1)>fabs(diff1))
			{
				min1=diff1;
				rxmin1=rx;		
			}
		}
		for (int i=1; i<=count-1; i++)
			exsum=exsum+pow(rx,i);							//building the geometric series
		diff2=L-exsum*deltax_start;							//CONSTRAINT 2: all cell-sizes SHOULD sum up to domain length	
		{										//store the minimum value of "diff2"
			if (fabs(min2)>fabs(diff2))
			{							
				min2=diff2;
				rxmin2=rx;
			}
		}
	if (rx<0.5 || rx>2.0)									//condition to prevent excessively stretched/compressed solutions
		break;
	else
		rx=rx+drx;
	}while((fabs(diff1)>tolerance) || (fabs(diff2)>tolerance)); 				//both constraints need to be simultaneously fulfilled	

	printf("count=%d \t diff1=%f \t rx1=%f \t diff2=%f \t rx2=%f\n",count,min1,rxmin1,min2,rxmin2);

if (rx>0.5 && rx<2.0)										//constraints fulfilled with acceptable stretching ratio
	{
		nx=count;
		rx=rx-((deltax_end-deltax_start)/fabs(deltax_end-deltax_start))*(mag_drx);	//to correct for the final increment/decrement in "rx"
		break;
	}											//solution is achieved only if rxmin1 = rxmin2
}
printf("\nThe requisite number of cells is %d and they need to be stretched by a factor of %f\n\n",nx,rx);
}
