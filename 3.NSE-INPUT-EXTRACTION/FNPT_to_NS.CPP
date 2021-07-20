//ROUTINE TO INTERPOLATE MOVING-GRID FNPT DATA ONTO A FIXED EULERIAN MESH FOR INPUT TO THE NAVIER-STOKES SOLVER 

#include <stdio.h>
#include <math.h>

main()
{
    FILE *f1;                   								//declaration of the file pointer for reading Output_PHIT.dat
    FILE *f2;                   								//declaration of the file pointer for mesh writing
    FILE *fU;                   								//declaration of the file pointer for writing to the Staggered U-mesh
    FILE *fW;                   								//declaration of the file pointer for writing to the Staggered W-mesh
    FILE *fp;                   								//declaration of the file pointer for writing to the Collocated p,f-mesh

    double d=1.0;               								//still water depth

    //input information for the NS domain
    double L=20.0;              								//length of the NS domain
    double H=2.0;               								//height of the NS domain
    int NX=500;                 								//total number of NS cells along x (has to be consistent with Global NS mesh size)
    int NZ=100;                 								//total number of NS cells along z (has to be consistent with Global NS mesh size)
    double dx=L/NX;             								//NS cell size along x
    double dz=H/NZ;             								//NS cell size along z

    //FNPT sub-domain extracted from IITM-FNPT2D - see "Output_PHIT.dat" and "FEMINPUT.dat"
    int NX_FNPT=10;             								//number of vertical FNPT layers along the x-direction [GRID IS FIXED]
    int NZ_FNPT=21;             								//number of horizontal FNPT layers along the z-direction [GRID MOVES VERTICALLY]

    //input information for the NS sub-domain where FNPT-data would be interpolated		/*...HARD-CODED - DO NOT CHANGE...*/
    int NX_NS=5;                								//number of cells along the length of the NS sub-domain
    int NZ_NS=NZ;               								//number of cells along the height of the NS sub-domain

    //number of time-steps executed in FNPT
    //N=(total no. of lines in "Output_PHIT.dat")/(NX_FNPT*NZ_FNPT)
    int N=5513;     										//verify against "FEMINPUT.dat" and "Output_PHIT.dat"

    //generation of the staggered U-mesh [back-staggered by 0.5*dx in the x-direction]
    double x_U[NX_NS+10]={0.0},z[NZ_NS+10]={0.0};
    x_U[1]=0.0;
    z[1]=0.5*dz;
    for (int i=2; i<=NX_NS; i++)								/***********************************************************************/
        x_U[i]=x_U[i-1]+dx;									/*...Change this routine if z-direction meshing is non-uniform in NS...*/
    for (int k=2; k<=NZ_NS; k++)								/***********************************************************************/
        z[k]=z[k-1]+dz;
    //writing the U-mesh to file (can be opened in TecPLOT)
    f2=fopen("Umesh.dat","w");
    fprintf (f2,"\nVARIABLES = ""x"", ""y""\n");
    fprintf (f2,"ZONE I=%d, J=%d, ZONETYPE=ORDERED, DATAPACKING=POINT\n",NZ_NS,NX_NS);
    for (int i=1; i<=NX_NS; i++)
        for (int k=1; k<=NZ_NS; k++)
            fprintf(f2,"%f\t%f\n",x_U[i],z[k]);
    fclose(f2);

    //generation of the collocated p,f-mesh
    double x[NX_NS+10]={0.0};
    x[1]=0.5*dx;
    z[1]=0.5*dz;
    for (int i=2; i<=NX_NS; i++)								/***********************************************************************/
        x[i]=x[i-1]+dx;										/*...Change this routine if z-direction meshing is non-uniform in NS...*/
    for (int k=2; k<=NZ_NS; k++)								/***********************************************************************/
        z[k]=z[k-1]+dz;
    //writing the p,f-mesh to file (can be opened in TecPLOT)
    f1=fopen("pmesh.dat","w");
    fprintf (f1,"\nVARIABLES = ""x"", ""y""\n");
    fprintf (f1,"ZONE I=%d, J=%d, ZONETYPE=ORDERED, DATAPACKING=POINT\n",NZ_NS,NX_NS);
    for (int i=1; i<=NX_NS; i++)
        for (int k=1; k<=NZ_NS; k++)
            fprintf(f1,"%f\t%f\n",x[i],z[k]);
    fclose(f1);

    //generation of the staggered W-mesh [back-staggered by 0.5*dz in the z-direction]
    double z_W[NZ_NS+10]={0.0};
    x[1]=0.5*dx;
    z_W[1]=0.0;
    for (int i=2; i<=NX_NS; i++)								/***********************************************************************/
        x[i]=x[i-1]+dx;										/*...Change this routine if z-direction meshing is non-uniform in NS...*/
    for (int k=2; k<=NZ_NS; k++)								/***********************************************************************/
        z_W[k]=z_W[k-1]+dz;
    //writing the W-mesh to file (can be opened in TecPLOT)
    f1=fopen("Wmesh.dat","w");
    fprintf (f1,"\nVARIABLES = ""x"", ""y""\n");
    fprintf (f1,"ZONE I=%d, J=%d, ZONETYPE=ORDERED, DATAPACKING=POINT\n",NZ_NS,NX_NS);
    for (int i=1; i<=NX_NS; i++)
        for (int k=1; k<=NZ_NS; k++)
            fprintf(f1,"%f\t%f\n",x[i],z_W[k]);
    fclose(f1);

    /*******************************************************************************************/	/**********************************/
    /*...Order of variables in Output_PHIT.dat : Time , X , Z , pressure , vx , vz , ax , az...*/	/*...HARD-CODED - DO NOT CHANGE...*/
    /*******************************************************************************************/	/**********************************/
    double FNPT_z[25]={0.0},FNPT_p[25]={0.0},FNPT_U[25]={0.0},FNPT_W[25]={0.0};
    f1=fopen("Output_ PHIT.dat","r");    							/*...open the FNPT output file in read mode...*/                    /*...KEEP IT OPEN...*/
    fU=fopen("Interp_Ustag.dat","w");   							/*...open a file to write U-velocity data over the U-grid...*/      /*...KEEP IT OPEN...*/
    fW=fopen("Interp_Wstag.dat","w");   							/*...open a file to write W-velocity data over the W-grid...*/      /*...KEEP IT OPEN...*/
    fp=fopen("Interp_pf.dat","w");      							/*...open a file to write p,f data over the collocated grid...*/    /*...KEEP IT OPEN...*/
    for (int count=1; count<=N; count++)    							//loop for all time steps [no time interpolation involved at this stage]
    {
        for (int i=1; i<=NX_FNPT/2; i++)    							//only half of the NX_FNPT nodes need to be considered for the U or pfW grids
        {
            /**************************************************************/
            /*...every (2*N-1) "FNPT_NZ lines" to be read to the U-grid...*/
            /**************************************************************/
            for (int k=1; k<=NZ_FNPT; k++)
            {
                fscanf(f1,"%*f\t%*f\t%lf\t%*f\t%lf\t%*f\t%*f\t%*f\n",&FNPT_z[k],&FNPT_U[k]); 	//only read z,p,vx,vz (skip the rest)
                FNPT_z[k]=FNPT_z[k]+d;  							//shifting the SWL to z=d
            }
                //*...interpolating FNPT information to the staggered U-grid...*//
                //*...interpolation is ONLY vertical...*//
                for (int l=1; l<=NZ_NS; l++)    						//vertical traverse over the staggered U-grid
                {
                    if (z[l]<FNPT_z[NZ_FNPT])   						//either a pure-cell in water or a mixed cell [U-cell center lies below the FNPT node]
                    {
                        for (int k=1; k<=NZ_FNPT; k++)
                        {
                        if (FNPT_z[k]>=z[l])
                        {
                            fprintf(fU,"%f\t%f\t%f\n",x_U[i],z[l],(((FNPT_U[k-1]*(FNPT_z[k]-z[l]))+(FNPT_U[k]*(z[l]-FNPT_z[k-1])))/(FNPT_z[k]-FNPT_z[k-1])));
                            break;
                        }
                        }
                    }
                    else                        						//either a pure-cell in air or a mixed cell [U-cell center lies above the FNPT node]
                            fprintf(fU,"%f\t%f\t%f\n",x_U[i],z[l],0.0);     			//surface velocity directly assigned
                }

            /*************************************************************************/
            /*...every (2*N) "FNPT_NZ lines" to be read to the W-grid and p,f-grid...*/
            /*************************************************************************/
            for (int k=1; k<=NZ_FNPT; k++)
            {
                fscanf(f1,"%*f\t%*f\t%lf\t%lf\t%*f\t%lf\t%*f\t%*f\n",&FNPT_z[k],&FNPT_p[k],&FNPT_W[k]); //only read z,p,vx,vz (skip rest)
                FNPT_z[k]=FNPT_z[k]+d;  							//shifting the SWL to z=d
            }
                //*...interpolating FNPT information to the staggered W-grid...*//
                //*...interpolation is ONLY vertical...*//
                for (int l=1; l<=NZ_NS; l++)    						//vertical traverse over the staggered W-grid
                {
                    if (z_W[l]<FNPT_z[NZ_FNPT]) 						//either a pure-cell in water or a mixed cell [W-cell center lies below the FNPT node]
                    {
                        for (int k=1; k<=NZ_FNPT; k++)
			{
                        if (FNPT_z[k]>=z_W[l])
                        {
			    if (FNPT_z[k]-FNPT_z[k-1]==0)
                            fprintf(fW,"%f\t%f\t%f\n",x[i],z_W[l],0.0);
			    else
                            fprintf(fW,"%f\t%f\t%f\n",x[i],z_W[l],(((FNPT_W[k-1]*(FNPT_z[k]-z_W[l]))+(FNPT_W[k]*(z_W[l]-FNPT_z[k-1])))/(FNPT_z[k]-FNPT_z[k-1])));
                            break;
                        }
			}
                    }
                    else                        						//either a pure-cell in air or a mixed cell [W-cell center lies above the FNPT node]
                            fprintf(fW,"%f\t%f\t%f\n",x[i],z_W[l],0.0);
                }
                //*...interpolating FNPT information to the collocated p,f-grid...*//
                //*...interpolation is ONLY vertical...*//
                for (int l=1; l<=NZ_NS; l++)    						//vertical traverse over the collocated p,f-grid
                {
                    if ((z[l]+0.5*dz)<=FNPT_z[NZ_FNPT])                                     	//pure cell lying entirely in water [p=rho*g*z and f=1]
                    {
                        for (int k=1; k<=NZ_FNPT; k++)
			{
                        if (z[l]>=FNPT_z[k])
                        {
                            fprintf(fp,"%f\t%f\t%f\t%f\n",x[i],z[l],1000.0*(((FNPT_p[k]*(FNPT_z[k+1]-z[l]))+(FNPT_p[k+1]*(z[l]-FNPT_z[k])))/(FNPT_z[k+1]-FNPT_z[k])),1.0);
                            break;
                        }
			}
                    }
                    else if (fabs(z[l]-FNPT_z[NZ_FNPT])<0.5*dz && FNPT_z[NZ_FNPT]>z[l])     	//mixed cell with the FNPT node lying above the p,f-cell center
                    {
                            int k=NZ_FNPT-1;
                            fprintf(fp,"%f\t%f\t%f\t%f\n",x[i],z[l],1000.0*(((FNPT_p[k]*(FNPT_z[k+1]-z[l]))+(FNPT_p[k+1]*(z[l]-FNPT_z[k])))/(FNPT_z[k+1]-FNPT_z[k])),(((FNPT_z[NZ_FNPT]-z[l])+0.5*dz)/dz));
                    }
                    else if (fabs(z[l]-FNPT_z[NZ_FNPT])<0.5*dz && FNPT_z[NZ_FNPT]<z[l])     	//mixed cell with the FNPT node lying below the p,f-cell center
                            fprintf(fp,"%f\t%f\t%f\t%f\n",x[i],z[l],0.0,(((FNPT_z[NZ_FNPT]-z[l])+0.5*dz)/dz));
                    else                                                                    	//pure cell lying entirely in air [p=0.0 and f=0.0]
                            fprintf(fp,"%f\t%f\t%f\t%f\n",x[i],z[l],0.0,0.0);
                }
        }//end i-loop
    }//end time-loop
    fclose(f1);
    fclose(fU);
    fclose(fW);
    fclose(fp);
}
