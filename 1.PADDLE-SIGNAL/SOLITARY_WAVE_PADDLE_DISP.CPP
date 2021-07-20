/**************************************************************************************************************************/
/*THIS CODE CALCULATES THE TOPOLOGY AS WELL AS PADDLE DISPLACEMENT CHARACTERISTICS REQUIRED FOR GENERATING A SOLITARY WAVE*/
/***********...GORING'S THEORY IS USED WHERE THE SOLITON HAS A HEIGHT "H" PROPAGATING IN A WATER DEPTH OF "d"...***********/
/**************************************************************************************************************************/

#include <stdio.h>
#include <math.h>

main()
{
double H=0.6;                               			//wave-height (m)
double d=1.0;                               			//water-depth (m)
double g=9.81;                              			//acceleration due to gravity
double C=sqrt(g*(H+d));                     			//phase velocity (m/s)
double k=sqrt((3.0*H)/(4.0*d*d*d));         			//wavenumber "k" (m^(-1))
double S=2.0*H/(k*d);                       			//required wave-maker stroke (m)
double T=(2.0*(3.8+(H/d)))/(k*C);           			//duration of soliton generation (s)
double T_init=2.0;                          			//initial period of inactivity (s)
double T_final=20.0;                       		 	//final period of inactivity (s)
double dt=0.01;                             			//step-size or sampling-rate

FILE *f1;                                   			//declaration of file-pointer
f1=fopen("wavePaddle.dat","a");             			//opening in a datafile in "append" mode

/*generation of time axis for the initial period of inactivity*/
double t_init[50000]={0.0};                 			//declaration of time array
int N_init=T_init/dt;                      			//number of time-steps required
t_init[0]=0.0;                              			//initial time
for (int i=1; i<=N_init; i++)
    t_init[i]=t_init[i-1]+dt;
for (int i=0; i<=N_init-1; i++)
    fprintf(f1,"%f\t%f\n",t_init[i],0.0);

/*generation of the time axis for period of soliton generation*/
double t[50000]={0.0};                      			//declaration of the time array
int N=T/dt;                                 			//number of time-steps required
t[0]=0.0;                                   			//initial time
for (int i=1; i<=N; i++)
    t[i]=t[i-1]+dt;                         			//generation of the time-vector

/*iterative calculation of the displacement*/
/*displacement formula modified so that WM only moves forwards in time*/
double X=0.0;                               			//initialization of wave-paddle displacement
double Xold=0.0;                            			//guessed value of paddle displacement
for (int i=0; i<=N; i++)                    			//loop over the entire time-vector
{
    do                                    			//simple RHS-LHS based iterative method
    {
        Xold=X;
        X=(H/(k*d))*(1.0+tanh(k*((C*(t[i]-0.5*T))-(Xold-0.5*S))));
        printf("Xold=%f\tX=%f\n",Xold,X);
    }while(fabs(X-Xold)>0.000000000001);

    do                                      			//Newton-Raphson method
    {
        Xold=X;
        X=Xold-(Xold-((H/(k*d))*(1.0+tanh(k*((C*(t[i]-0.5*T))-(Xold-0.5*S))))))/(1.0+(H/d)/(pow((cosh((k*C*(t[i]-0.5*T))-(k*(Xold-0.5*S)))),2.0)));
        printf("Xold=%f\tX=%f\n",Xold,X);
    }while(fabs(X-Xold)>0.000000000001);

    fprintf(f1,"%f\t%f\n",t[i]+T_init,X);
}

/*generation of time axis for the final period of inactivity*/
double t_final[50000]={0.0};                			//declaration of time array
int N_final=T_final/dt;
t_final[0]=T_init+t[N]+dt;
for (int i=1; i<=N_final; i++)
    t_final[i]=t_final[i-1]+dt;
for (int i=0; i<=N_final-1; i++)
    fprintf(f1,"%f\t%f\n",t_final[i],X);

fclose(f1);
}
