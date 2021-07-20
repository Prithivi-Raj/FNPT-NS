#define NX 500											/***************************/
#define NY 50											/*...Global NS mesh size...*/
#define NZ 100											/***************************/

/****************************************************/
/*...HARD-CODED FUNCTIONS CALLED BY THE MAIN CODE...*/
/*.....DO NOT MAKE ANY CHANGES BELOW THIS POINT.....*/
/****************************************************/

/*****************************************************************************/
/*....................NSE SOLVER-specific declarations.......................*/
/*****************************************************************************/
int i,j,k;
int IMG,sch,count,gcount,N_water;
double L,B,H,rx,ry,rz;
double exsum;
double ue,uw,un,us,ut,ub,ve,vw,vn,vs,vt,vb,we,ww,wn,ws,wt,wb;
double uea,uwa,una,usa,uta,uba,vea,vwa,vna,vsa,vta,vba,wea,wwa,wna,wsa,wta,wba;
double uadvel,phiuU,phiuUU,phiuUUU,phiuD,phiuDD;
double vadvel,phivU,phivUU,phivUUU,phivD,phivDD;
double wadvel,phiwU,phiwUU,phiwUUU,phiwD,phiwDD;
double VOLU,VOLUU,VOLUUU,VOLD,VOLDD;
double lx[NX+10],ly[NY+10],lz[NZ+10],cx[NX+10],cy[NY+10],cz[NZ+10];
double u[NX+10][NY+10][NZ+10],v[NX+10][NY+10][NZ+10],w[NX+10][NY+10][NZ+10],utild[NX+10][NY+10][NZ+10],vtild[NX+10][NY+10][NZ+10],wtild[NX+10][NY+10][NZ+10],ustar[NX+10][NY+10][NZ+10],vstar[NX+10][NY+10][NZ+10],wstar[NX+10][NY+10][NZ+10];
double p[NX+10][NY+10][NZ+10],pcorr[NX+10][NY+10][NZ+10];
double must[NX+10][NY+10][NZ+10],rhost[NX+10][NY+10][NZ+10];
double Uo,rho_liquid,rho_gas,rho_solid,mu_liquid,mu_gas,mu_solid,g;
double rhoe,rhow,rhon,rhos,rhot,rhob,mue,muw,mun,mus,mut,mub;
double dxe,dxw,dyn,dys,dzt,dzb;
double mass[NX+10][NY+10][NZ+10],gsum,apd,res[NX+10][NY+10][NZ+10],corr,ap,gsq,gRMS;
double usum,vsum,wsum,OLD_U,OLD_V,OLD_W,udiff,usq,uRMS,vdiff,vsq,vRMS,wdiff,wsq,wRMS;
double Unew[21][21][NZ+10]={0.0},Wnew[21][21][NZ+10]={0.0},Pnew[21][21][NZ+10]={0.0},Fnew[21][21][NZ+10]={0.0};//for storing the new FNPT field
double Uold[21][21][NZ+10]={0.0},Wold[21][21][NZ+10]={0.0},Pold[21][21][NZ+10]={0.0},Fold[21][21][NZ+10]={0.0};//for storing the old FNPT field
double tme,OLD_tme_FNPT,NEW_tme_FNPT,dt,dt_FNPT;
double FXAD1,FXAD2,FXAD,FYAD1,FYAD2,FYAD,FZAD1,FZAD2,FZAD,FXDF,FYDF,FZDF,FZSRC_GR;
char string[180];
double uplot[NX+10][NY+10][NZ+10],vplot[NX+10][NY+10][NZ+10],wplot[NX+10][NY+10][NZ+10];
double Pnew_fnpt=0.0,Wnew_fnpt=0.0,Pold_fnpt=0.0,Wold_fnpt=0.0,p_fnpt=0.0,w_fnpt=0.0;		//FNPT and NS p,w fields to be compared at the first cell
double water_depth,EOPC_ERROR;

/*****************************************************************************/
/*......................CICSAM-specific declarations.........................*/
/*****************************************************************************/
int swp;
double f[NX+10][NY+10][NZ+10]={0.0},fnew[NX+10][NY+10][NZ+10]={0.0};				//general volume fractions
double c[NX+10][NY+10][NZ+10]={0.0};								//general color function
int CICSAM_RUN;											//FFD treatment of three phases
double liquid_VOF[NX+10][NY+10][NZ+10]={0.0},solid_VOF[NX+10][NY+10][NZ+10]={0.0},initsolid[NX+10][NY+10][NZ+10]={0.0};//FFD treatment of three phases
double epsilon=0.0000000000001;
double pi=3.141592653589793238462643383279502884197;
int forlorn,N;
int a1,b1,c1,usc=0,osc=0;
double min,max,temp,Dmin,Dmax;
double vel,flux,len,s1,s2;
double fe,fw,fn,fs,ft,fb;
double initVOL,finVOL;
double face_COUR,f_DONOR,f_UPWIND,f_ACCEPT;
double ftild_DONOR,ftild_face_CBC,ftild_face_UQ,ftild_face,theta_f,gamma_f,beta_f,norx,nory,norz;
double f_U,f_UU,f_D;

/***********************************************/
/*...BASIC CONSTITUTION OF THE CICSAM SCHEME...*/
/***********************************************/
void stages_CICSAM()
{
            //CICSAM [First Component] ----> HYPER-C [compressive scheme]
            {
                if (ftild_DONOR<(0.0-epsilon) || ftild_DONOR>(1.0+epsilon))
                    ftild_face_CBC=ftild_DONOR;

                else if (ftild_DONOR>=(0.0-epsilon) && ftild_DONOR<=(1.0+epsilon))
                {
                    if (face_COUR>0.0+epsilon)
                    {
                    if ((ftild_DONOR/face_COUR)<1.0)
                        ftild_face_CBC=ftild_DONOR/face_COUR;
                    else
                        ftild_face_CBC=1.0;
                    }
                    else
                        ftild_face_CBC=ftild_DONOR;						//to ensure face_VOF * vel = 0.0
                }
                else
                    ftild_face_CBC=ftild_DONOR;							//fall-back to UPWINDING
            }

            //CICSAM [Second Component] ----> Ultimate-QUICKEST [diffusive scheme]
            {
                if (ftild_DONOR<(0.0-epsilon) || ftild_DONOR>(1.0+epsilon))
                    ftild_face_UQ=ftild_DONOR;

                else if (ftild_DONOR>=(0.0-epsilon) && ftild_DONOR<=(1.0+epsilon))
                {
                    if (face_COUR>0.0+epsilon)
                    {
                        if (((8.0*face_COUR*ftild_DONOR+(1.0-face_COUR)*(6.0*ftild_DONOR+3.0))/(8.0))<ftild_face_CBC)
                            ftild_face_UQ=((8.0*face_COUR*ftild_DONOR+(1.0-face_COUR)*(6.0*ftild_DONOR+3.0))/(8.0));
                        else
                            ftild_face_UQ=ftild_face_CBC;
                    }
                    else
                        ftild_face_UQ=ftild_DONOR;
                }
                else
                    ftild_face_UQ=ftild_DONOR;							//fall-back to UPWINDING
            }

            //blending between HYPER-C and Ultimate-QUICKEST
            {
                if ((0.5*(1.0+cos(2.0*theta_f)))<1.0)
                    gamma_f=(0.5*(1.0+cos(2.0*theta_f)));
                else
                    gamma_f=1.0;
            }

            ftild_face=gamma_f*ftild_face_CBC+(1.0-gamma_f)*ftild_face_UQ;
            if (ftild_DONOR>=1.0-epsilon && ftild_DONOR<=1.0+epsilon)
                beta_f=0.0;									//switch to UPWINDING if f_DONOR = f_ACCEPT
            else if (fabs(f_DONOR-f_ACCEPT)<=epsilon && fabs(f_DONOR-f_UPWIND)<=epsilon)
                beta_f=0.0;									//switch to UPWINDING if f_DONOR = f_ACCEPT = f_UPWIND
            else if (fabs(f_ACCEPT-f_UPWIND)<=epsilon)
                beta_f=0.0;									//switch to UPWINDING if f_ACCEPT = f_UPWIND
            else
                beta_f=(ftild_face-ftild_DONOR)/(1.0-ftild_DONOR);
}

/********************************************/
/*...IMPLEMENTING CICSAM TO THE VOF FIELD...*/
/********************************************/
void f_CICSAM()
{
    if (swp==0)											//BEGIN X-SWEEP.................................................
    {
        if (ue>0)										//fluxing is OUTWARD (ue is +ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j][k];
            f_ACCEPT=f[i+1][j][k];
            f_UPWIND=f[i-1][j][k];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j,k]
            double ue_D,uw_D;
            face_COUR=0.0;
            {
            ue_D=(lx[i+1]*u[i][j][k]+lx[i]*u[i+1][j][k])/(lx[i+1]+lx[i]);
            uw_D=(lx[i]*u[i-1][j][k]+lx[i-1]*u[i][j][k])/(lx[i]+lx[i-1]);
                face_COUR=fmax(ue_D*dt/lx[i],0.0)+fmax(-uw_D*dt/lx[i],0.0);
            }

            //calculation of the interface normal (Parker and Youngs' method)
            double E_nt=((f[i+1][j+1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_nt=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_nt=((f[i-1][j+1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_st=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_st=((f[i-1][j][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_st=((f[i+1][j][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_nb=((f[i+1][j+1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_nb=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_nb=((f[i-1][j+1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_sb=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_sb=((f[i-1][j][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_sb=((f[i+1][j][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_et=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_wt=((f[i][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_wt=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_wt=((f[i][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_et=((f[i+1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_eb=((f[i+1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_eb=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_wb=((f[i][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_wb=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_wb=((f[i][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_eb=((f[i+1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_ne=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double T_nw=((f[i][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_nw=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_sw=((f[i][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_sw=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double T_se=((f[i+1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_se=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_ne=((f[i+1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_nw=((f[i][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_sw=((f[i][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_se=((f[i+1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(norx/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along +x

            stages_CICSAM();

            fe=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }

        else											//fluxing is INWARD (ue is -ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i+1][j][k];
            f_ACCEPT=f[i][j][k];
            f_UPWIND=f[i+2][j][k];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i+1,j,k]
            double ue_D,uw_D;
            face_COUR=0.0;
            {
            ue_D=(lx[i+2]*u[i+1][j][k]+lx[i+1]*u[i+2][j][k])/(lx[i+2]+lx[i+1]);
            uw_D=(lx[i+1]*u[i][j][k]+lx[i]*u[i+1][j][k])/(lx[i+1]+lx[i]);
                face_COUR=fmax(ue_D*dt/lx[i],0.0)+fmax(-uw_D*dt/lx[i],0.0);
            }

            /*...calculation of the interface normal (Parker and Youngs' method) [stencil shifts one cell along +i]...*/
            double E_nt=((f[i+2][j+1][k+1]*ly[j]*lz[k])+(f[i+2][j][k]*ly[j+1]*lz[k+1])+(f[i+2][j][k+1]*ly[j+1]*lz[k])+(f[i+2][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_nt=((f[i+1][j+1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_nt=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_st=((f[i+1][j][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_st=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_st=((f[i+2][j][k+1]*ly[j-1]*lz[k])+(f[i+2][j-1][k]*ly[j]*lz[k+1])+(f[i+2][j-1][k+1]*ly[j]*lz[k])+(f[i+2][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_nb=((f[i+2][j+1][k]*ly[j]*lz[k-1])+(f[i+2][j][k-1]*ly[j+1]*lz[k])+(f[i+2][j][k]*ly[j+1]*lz[k-1])+(f[i+2][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_nb=((f[i+1][j+1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_nb=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_sb=((f[i+1][j][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_sb=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_sb=((f[i+2][j][k]*ly[j-1]*lz[k-1])+(f[i+2][j-1][k-1]*ly[j]*lz[k])+(f[i+2][j-1][k]*ly[j]*lz[k-1])+(f[i+2][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+2]+lx[i+1]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_3=(P_st-W_st)/(0.5*(lx[i+1]+lx[i]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+2]+lx[i+1]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+2]+lx[i+1]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i+1]+lx[i]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+2]+lx[i+1]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+2][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i+2]*lz[k+1])+(f[i+1][j+1][k+1]*lx[i+2]*lz[k])+(f[i+2][j+1][k]*lx[i+1]*lz[k+1]))/((lx[i+1]+lx[i+2])*(lz[k]+lz[k+1]));
            double P_et=((f[i+2][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i+2]*lz[k+1])+(f[i+1][j][k+1]*lx[i+2]*lz[k])+(f[i+2][j][k]*lx[i+1]*lz[k+1]))/((lx[i+1]+lx[i+2])*(lz[k]+lz[k+1]));
            double N_wt=((f[i+1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_wt=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double S_wt=((f[i+1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double S_et=((f[i+2][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i+2]*lz[k+1])+(f[i+1][j-1][k+1]*lx[i+2]*lz[k])+(f[i+2][j-1][k]*lx[i+1]*lz[k+1]))/((lx[i+1]+lx[i+2])*(lz[k]+lz[k+1]));
            double N_eb=((f[i+2][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i+2]*lz[k])+(f[i+1][j+1][k]*lx[i+2]*lz[k-1])+(f[i+2][j+1][k-1]*lx[i+1]*lz[k]))/((lx[i+1]+lx[i+2])*(lz[k-1]+lz[k]));
            double P_eb=((f[i+2][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i+2]*lz[k])+(f[i+1][j][k]*lx[i+2]*lz[k-1])+(f[i+2][j][k-1]*lx[i+1]*lz[k]))/((lx[i+1]+lx[i+2])*(lz[k-1]+lz[k]));
            double N_wb=((f[i+1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_wb=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double S_wb=((f[i+1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double S_eb=((f[i+2][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i+2]*lz[k])+(f[i+1][j-1][k]*lx[i+2]*lz[k-1])+(f[i+2][j-1][k-1]*lx[i+1]*lz[k]))/((lx[i+1]+lx[i+2])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+2][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i+2]*ly[j+1])+(f[i+1][j+1][k+1]*lx[i+2]*ly[j])+(f[i+2][j][k+1]*lx[i+1]*ly[j+1]))/((lx[i+1]+lx[i+2])*(ly[j]+ly[j+1]));
            double P_ne=((f[i+2][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i+2]*ly[j+1])+(f[i+1][j+1][k]*lx[i+2]*ly[j])+(f[i+2][j][k]*lx[i+1]*ly[j+1]))/((lx[i+1]+lx[i+2])*(ly[j]+ly[j+1]));
            double T_nw=((f[i+1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_nw=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double T_sw=((f[i+1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_sw=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double T_se=((f[i+2][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i+2]*ly[j])+(f[i+1][j][k+1]*lx[i+2]*ly[j-1])+(f[i+2][j-1][k+1]*lx[i+1]*ly[j]))/((lx[i+1]+lx[i+2])*(ly[j-1]+ly[j]));
            double P_se=((f[i+2][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i+2]*ly[j])+(f[i+1][j][k]*lx[i+2]*ly[j-1])+(f[i+2][j-1][k]*lx[i+1]*ly[j]))/((lx[i+1]+lx[i+2])*(ly[j-1]+ly[j]));
            double B_ne=((f[i+2][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i+2]*ly[j+1])+(f[i+1][j+1][k-1]*lx[i+2]*ly[j])+(f[i+2][j][k-1]*lx[i+1]*ly[j+1]))/((lx[i+1]+lx[i+2])*(ly[j]+ly[j+1]));
            double B_nw=((f[i+1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_sw=((f[i+1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_se=((f[i+2][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i+2]*ly[j])+(f[i+1][j][k-1]*lx[i+2]*ly[j-1])+(f[i+2][j-1][k-1]*lx[i+1]*ly[j]))/((lx[i+1]+lx[i+2])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(-norx/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along -x

            stages_CICSAM();

            fe=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }

        if (uw<0)										//fluxing is OUTWARD (uw is -ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j][k];
            f_ACCEPT=f[i-1][j][k];
            f_UPWIND=f[i+1][j][k];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j,k]
            double ue_D,uw_D;
            face_COUR=0.0;
            {
            ue_D=(lx[i+1]*u[i][j][k]+lx[i]*u[i+1][j][k])/(lx[i+1]+lx[i]);
            uw_D=(lx[i]*u[i-1][j][k]+lx[i-1]*u[i][j][k])/(lx[i]+lx[i-1]);
                face_COUR=fmax(ue_D*dt/lx[i],0.0)+fmax(-uw_D*dt/lx[i],0.0);
            }

            //calculation of the interface normal (Parker and Youngs' method)
            double E_nt=((f[i+1][j+1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_nt=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_nt=((f[i-1][j+1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_st=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_st=((f[i-1][j][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_st=((f[i+1][j][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_nb=((f[i+1][j+1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_nb=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_nb=((f[i-1][j+1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_sb=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_sb=((f[i-1][j][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_sb=((f[i+1][j][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_et=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_wt=((f[i][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_wt=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_wt=((f[i][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_et=((f[i+1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_eb=((f[i+1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_eb=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_wb=((f[i][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_wb=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_wb=((f[i][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_eb=((f[i+1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_ne=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double T_nw=((f[i][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_nw=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_sw=((f[i][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_sw=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double T_se=((f[i+1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_se=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_ne=((f[i+1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_nw=((f[i][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_sw=((f[i][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_se=((f[i+1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(-norx/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along -x

            stages_CICSAM();

            fw=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }

        else											//fluxing is INWARD (uw is +ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i-1][j][k];
            f_ACCEPT=f[i][j][k];
            f_UPWIND=f[i-2][j][k];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i-1,j,k]
            double ue_D,uw_D;
            face_COUR=0.0;
            {
            ue_D=(lx[i]*u[i-1][j][k]+lx[i-1]*u[i][j][k])/(lx[i]+lx[i-1]);
            uw_D=(lx[i-1]*u[i-2][j][k]+lx[i-2]*u[i-1][j][k])/(lx[i-1]+lx[i-2]);
                face_COUR=fmax(ue_D*dt/lx[i],0.0)+fmax(-uw_D*dt/lx[i],0.0);
            }

            /*...calculation of the interface normal (Parker and Youngs' method) [stencil shifts one cell along -i]...*/
            double E_nt=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_nt=((f[i-1][j+1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_nt=((f[i-2][j+1][k+1]*ly[j]*lz[k])+(f[i-2][j][k]*ly[j+1]*lz[k+1])+(f[i-2][j][k+1]*ly[j+1]*lz[k])+(f[i-2][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_st=((f[i-1][j][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_st=((f[i-2][j][k+1]*ly[j-1]*lz[k])+(f[i-2][j-1][k]*ly[j]*lz[k+1])+(f[i-2][j-1][k+1]*ly[j]*lz[k])+(f[i-2][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_st=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_nb=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_nb=((f[i-1][j+1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_nb=((f[i-2][j+1][k]*ly[j]*lz[k-1])+(f[i-2][j][k-1]*ly[j+1]*lz[k])+(f[i-2][j][k]*ly[j+1]*lz[k-1])+(f[i-2][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_sb=((f[i-1][j][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_sb=((f[i-2][j][k]*ly[j-1]*lz[k-1])+(f[i-2][j-1][k-1]*ly[j]*lz[k])+(f[i-2][j-1][k]*ly[j]*lz[k-1])+(f[i-2][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_sb=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i-1]+lx[i-2]));
            double term_3=(P_st-W_st)/(0.5*(lx[i-1]+lx[i-2]));
            double term_4=(E_st-P_st)/(0.5*(lx[i]+lx[i-1]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i-1]+lx[i-2]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i-1]+lx[i-2]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i]+lx[i-1]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_et=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double N_wt=((f[i-1][j+1][k+1]*lx[i-2]*lz[k])+(f[i-2][j+1][k]*lx[i-1]*lz[k+1])+(f[i-2][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i-2]*lz[k+1]))/((lx[i-2]+lx[i-1])*(lz[k]+lz[k+1]));
            double P_wt=((f[i-1][j][k+1]*lx[i-2]*lz[k])+(f[i-2][j][k]*lx[i-1]*lz[k+1])+(f[i-2][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i-2]*lz[k+1]))/((lx[i-2]+lx[i-1])*(lz[k]+lz[k+1]));
            double S_wt=((f[i-1][j-1][k+1]*lx[i-2]*lz[k])+(f[i-2][j-1][k]*lx[i-1]*lz[k+1])+(f[i-2][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i-2]*lz[k+1]))/((lx[i-2]+lx[i-1])*(lz[k]+lz[k+1]));
            double S_et=((f[i][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double N_eb=((f[i][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_eb=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double N_wb=((f[i-1][j+1][k]*lx[i-2]*lz[k-1])+(f[i-2][j+1][k-1]*lx[i-1]*lz[k])+(f[i-2][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i-2]*lz[k]))/((lx[i-2]+lx[i-1])*(lz[k-1]+lz[k]));
            double P_wb=((f[i-1][j][k]*lx[i-2]*lz[k-1])+(f[i-2][j][k-1]*lx[i-1]*lz[k])+(f[i-2][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i-2]*lz[k]))/((lx[i-2]+lx[i-1])*(lz[k-1]+lz[k]));
            double S_wb=((f[i-1][j-1][k]*lx[i-2]*lz[k-1])+(f[i-2][j-1][k-1]*lx[i-1]*lz[k])+(f[i-2][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i-2]*lz[k]))/((lx[i-2]+lx[i-1])*(lz[k-1]+lz[k]));
            double S_eb=((f[i][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_ne=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_nw=((f[i-1][j+1][k+1]*lx[i-2]*ly[j])+(f[i-2][j][k+1]*lx[i-1]*ly[j+1])+(f[i-2][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i-2]*ly[j+1]))/((lx[i-2]+lx[i-1])*(ly[j]+ly[j+1]));
            double P_nw=((f[i-1][j+1][k]*lx[i-2]*ly[j])+(f[i-2][j][k]*lx[i-1]*ly[j+1])+(f[i-2][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i-2]*ly[j+1]))/((lx[i-2]+lx[i-1])*(ly[j]+ly[j+1]));
            double T_sw=((f[i-1][j][k+1]*lx[i-2]*ly[j-1])+(f[i-2][j-1][k+1]*lx[i-1]*ly[j])+(f[i-2][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i-2]*ly[j]))/((lx[i-2]+lx[i-1])*(ly[j-1]+ly[j]));
            double P_sw=((f[i-1][j][k]*lx[i-2]*ly[j-1])+(f[i-2][j-1][k]*lx[i-1]*ly[j])+(f[i-2][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i-2]*ly[j]))/((lx[i-2]+lx[i-1])*(ly[j-1]+ly[j]));
            double T_se=((f[i][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_se=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_ne=((f[i][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_nw=((f[i-1][j+1][k-1]*lx[i-2]*ly[j])+(f[i-2][j][k-1]*lx[i-1]*ly[j+1])+(f[i-2][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i-2]*ly[j+1]))/((lx[i-2]+lx[i-1])*(ly[j]+ly[j+1]));
            double B_sw=((f[i-1][j][k-1]*lx[i-2]*ly[j-1])+(f[i-2][j-1][k-1]*lx[i-1]*ly[j])+(f[i-2][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i-2]*ly[j]))/((lx[i-2]+lx[i-1])*(ly[j-1]+ly[j]));
            double B_se=((f[i][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(norx/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along +x

            stages_CICSAM();

            fw=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }
    }												//END X-SWEEP...................................................

    else if (swp==1)										//BEGIN Y-SWEEP.................................................
    {
        if (vn>0)										//fluxing is OUTWARD (vn is +ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j][k];
            f_ACCEPT=f[i][j+1][k];
            f_UPWIND=f[i][j-1][k];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j,k]
            double vn_D,vs_D;
            face_COUR=0.0;
            {
            vn_D=(ly[j+1]*v[i][j][k]+ly[j]*v[i][j+1][k])/(ly[j+1]+ly[j]);
            vs_D=(ly[j]*v[i][j-1][k]+ly[j-1]*v[i][j][k])/(ly[j]+ly[j-1]);
                face_COUR=fmax(vn_D*dt/ly[j],0.0)+fmax(-vs_D*dt/ly[j],0.0);
            }

            //calculation of the interface normal (Parker and Youngs' method)
            double E_nt=((f[i+1][j+1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_nt=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_nt=((f[i-1][j+1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_st=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_st=((f[i-1][j][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_st=((f[i+1][j][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_nb=((f[i+1][j+1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_nb=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_nb=((f[i-1][j+1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_sb=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_sb=((f[i-1][j][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_sb=((f[i+1][j][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_et=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_wt=((f[i][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_wt=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_wt=((f[i][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_et=((f[i+1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_eb=((f[i+1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_eb=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_wb=((f[i][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_wb=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_wb=((f[i][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_eb=((f[i+1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_ne=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double T_nw=((f[i][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_nw=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_sw=((f[i][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_sw=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double T_se=((f[i+1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_se=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_ne=((f[i+1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_nw=((f[i][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_sw=((f[i][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_se=((f[i+1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(nory/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along +y

            stages_CICSAM();

            fn=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }

        else											//fluxing is INWARD (vn is -ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j+1][k];
            f_ACCEPT=f[i][j][k];
            f_UPWIND=f[i][j+2][k];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j+1,k]
            double vn_D,vs_D;
            face_COUR=0.0;
            {
            vn_D=(ly[j+2]*v[i][j+1][k]+ly[j+1]*v[i][j+2][k])/(ly[j+2]+ly[j+1]);
            vs_D=(ly[j+1]*v[i][j][k]+ly[j]*v[i][j+1][k])/(ly[j+1]+ly[j]);
                face_COUR=fmax(vn_D*dt/ly[j],0.0)+fmax(-vs_D*dt/ly[j],0.0);
            }

            /*...calculation of the interface normal (Parker and Youngs' method) [stencil shifts one cell along +j]...*/
            double E_nt=((f[i+1][j+2][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j+2]*lz[k+1])+(f[i+1][j+1][k+1]*ly[j+2]*lz[k])+(f[i+1][j+2][k]*ly[j+1]*lz[k+1]))/((ly[j+1]+ly[j+2])*(lz[k]+lz[k+1]));
            double P_nt=((f[i][j+2][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j+2]*lz[k+1])+(f[i][j+1][k+1]*ly[j+2]*lz[k])+(f[i][j+2][k]*ly[j+1]*lz[k+1]))/((ly[j+1]+ly[j+2])*(lz[k]+lz[k+1]));
            double W_nt=((f[i-1][j+2][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j+2]*lz[k+1])+(f[i-1][j+1][k+1]*ly[j+2]*lz[k])+(f[i-1][j+2][k]*ly[j+1]*lz[k+1]))/((ly[j+1]+ly[j+2])*(lz[k]+lz[k+1]));
            double P_st=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_st=((f[i-1][j+1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double E_st=((f[i+1][j+1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double E_nb=((f[i+1][j+2][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j+2]*lz[k])+(f[i+1][j+1][k]*ly[j+2]*lz[k-1])+(f[i+1][j+2][k-1]*ly[j+1]*lz[k]))/((ly[j+1]+ly[j+2])*(lz[k-1]+lz[k]));
            double P_nb=((f[i][j+2][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j+2]*lz[k])+(f[i][j+1][k]*ly[j+2]*lz[k-1])+(f[i][j+2][k-1]*ly[j+1]*lz[k]))/((ly[j+1]+ly[j+2])*(lz[k-1]+lz[k]));
            double W_nb=((f[i-1][j+2][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j+2]*lz[k])+(f[i-1][j+1][k]*ly[j+2]*lz[k-1])+(f[i-1][j+2][k-1]*ly[j+1]*lz[k]))/((ly[j+1]+ly[j+2])*(lz[k-1]+lz[k]));
            double P_sb=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_sb=((f[i-1][j+1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double E_sb=((f[i+1][j+1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j+2][k+1]*lx[i]*lz[k])+(f[i][j+2][k]*lx[i+1]*lz[k+1])+(f[i][j+2][k+1]*lx[i+1]*lz[k])+(f[i+1][j+2][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_et=((f[i+1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_wt=((f[i][j+2][k+1]*lx[i-1]*lz[k])+(f[i-1][j+2][k]*lx[i]*lz[k+1])+(f[i-1][j+2][k+1]*lx[i]*lz[k])+(f[i][j+2][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_wt=((f[i][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_wt=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_et=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_eb=((f[i+1][j+2][k]*lx[i]*lz[k-1])+(f[i][j+2][k-1]*lx[i+1]*lz[k])+(f[i][j+2][k]*lx[i+1]*lz[k-1])+(f[i+1][j+2][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_eb=((f[i+1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_wb=((f[i][j+2][k]*lx[i-1]*lz[k-1])+(f[i-1][j+2][k-1]*lx[i]*lz[k])+(f[i-1][j+2][k]*lx[i]*lz[k-1])+(f[i][j+2][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_wb=((f[i][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_wb=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_eb=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+2]+ly[j+1]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+2]+ly[j+1]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_12=(P_et-S_et)/(0.5*(ly[j+1]+ly[j]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+2]+ly[j+1]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+2]+ly[j+1]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j+1]+ly[j]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j+2][k+1]*lx[i]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j+2])+(f[i][j+2][k+1]*lx[i+1]*ly[j+1])+(f[i+1][j+1][k+1]*lx[i]*ly[j+2]))/((lx[i]+lx[i+1])*(ly[j+1]+ly[j+2]));
            double P_ne=((f[i+1][j+2][k]*lx[i]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j+2])+(f[i][j+2][k]*lx[i+1]*ly[j+1])+(f[i+1][j+1][k]*lx[i]*ly[j+2]))/((lx[i]+lx[i+1])*(ly[j+1]+ly[j+2]));
            double T_nw=((f[i][j+2][k+1]*lx[i-1]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j+2])+(f[i-1][j+2][k+1]*lx[i]*ly[j+1])+(f[i][j+1][k+1]*lx[i-1]*ly[j+2]))/((lx[i-1]+lx[i])*(ly[j+1]+ly[j+2]));
            double P_nw=((f[i][j+2][k]*lx[i-1]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j+2])+(f[i-1][j+2][k]*lx[i]*ly[j+1])+(f[i][j+1][k]*lx[i-1]*ly[j+2]))/((lx[i-1]+lx[i])*(ly[j+1]+ly[j+2]));
            double T_sw=((f[i][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_sw=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_se=((f[i+1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_se=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_ne=((f[i+1][j+2][k-1]*lx[i]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j+2])+(f[i][j+2][k-1]*lx[i+1]*ly[j+1])+(f[i+1][j+1][k-1]*lx[i]*ly[j+2]))/((lx[i]+lx[i+1])*(ly[j+1]+ly[j+2]));
            double B_nw=((f[i][j+2][k-1]*lx[i-1]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j+2])+(f[i-1][j+2][k-1]*lx[i]*ly[j+1])+(f[i][j+1][k-1]*lx[i-1]*ly[j+2]))/((lx[i-1]+lx[i])*(ly[j+1]+ly[j+2]));
            double B_sw=((f[i][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_se=((f[i+1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(-nory/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along -y

            stages_CICSAM();

            fn=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }

        if (vs<0)										//fluxing is OUTWARD (vs is -ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j][k];
            f_ACCEPT=f[i][j-1][k];
            f_UPWIND=f[i][j+1][k];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j,k]
            double vn_D,vs_D;
            face_COUR=0.0;
            {
            vn_D=(ly[j+1]*v[i][j][k]+ly[j]*v[i][j+1][k])/(ly[j+1]+ly[j]);
            vs_D=(ly[j]*v[i][j-1][k]+ly[j-1]*v[i][j][k])/(ly[j]+ly[j-1]);
                face_COUR=fmax(vn_D*dt/ly[j],0.0)+fmax(-vs_D*dt/ly[j],0.0);
            }

            //calculation of the interface normal (Parker and Youngs' method)
            double E_nt=((f[i+1][j+1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_nt=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_nt=((f[i-1][j+1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_st=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_st=((f[i-1][j][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_st=((f[i+1][j][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_nb=((f[i+1][j+1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_nb=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_nb=((f[i-1][j+1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_sb=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_sb=((f[i-1][j][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_sb=((f[i+1][j][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_et=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_wt=((f[i][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_wt=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_wt=((f[i][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_et=((f[i+1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_eb=((f[i+1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_eb=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_wb=((f[i][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_wb=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_wb=((f[i][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_eb=((f[i+1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_ne=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double T_nw=((f[i][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_nw=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_sw=((f[i][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_sw=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double T_se=((f[i+1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_se=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_ne=((f[i+1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_nw=((f[i][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_sw=((f[i][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_se=((f[i+1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(-nory/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along -y

            stages_CICSAM();

            fs=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }

        else											//fluxing is INWARD (vs is +ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j-1][k];
            f_ACCEPT=f[i][j][k];
            f_UPWIND=f[i][j-2][k];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j-1,k]
            double vn_D,vs_D;
            face_COUR=0.0;
            {
            vn_D=(ly[j]*v[i][j-1][k]+ly[j-1]*v[i][j][k])/(ly[j]+ly[j-1]);
            vs_D=(ly[j-1]*v[i][j-2][k]+ly[j-2]*v[i][j-1][k])/(ly[j-1]+ly[j-2]);
                face_COUR=fmax(vn_D*dt/ly[j],0.0)+fmax(-vs_D*dt/ly[j],0.0);
            }

            /*...calculation of the interface normal (Parker and Youngs' method) [stencil shifts one cell along -j]...*/
            double E_nt=((f[i+1][j][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double P_nt=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_nt=((f[i-1][j][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double P_st=((f[i][j-1][k+1]*ly[j-2]*lz[k])+(f[i][j-2][k]*ly[j-1]*lz[k+1])+(f[i][j-2][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j-2]*lz[k+1]))/((ly[j-2]+ly[j-1])*(lz[k]+lz[k+1]));
            double W_st=((f[i-1][j-1][k+1]*ly[j-2]*lz[k])+(f[i-1][j-2][k]*ly[j-1]*lz[k+1])+(f[i-1][j-2][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j-2]*lz[k+1]))/((ly[j-2]+ly[j-1])*(lz[k]+lz[k+1]));
            double E_st=((f[i+1][j-1][k+1]*ly[j-2]*lz[k])+(f[i+1][j-2][k]*ly[j-1]*lz[k+1])+(f[i+1][j-2][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j-2]*lz[k+1]))/((ly[j-2]+ly[j-1])*(lz[k]+lz[k+1]));
            double E_nb=((f[i+1][j][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double P_nb=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_nb=((f[i-1][j][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double P_sb=((f[i][j-1][k]*ly[j-2]*lz[k-1])+(f[i][j-2][k-1]*ly[j-1]*lz[k])+(f[i][j-2][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j-2]*lz[k]))/((ly[j-2]+ly[j-1])*(lz[k-1]+lz[k]));
            double W_sb=((f[i-1][j-1][k]*ly[j-2]*lz[k-1])+(f[i-1][j-2][k-1]*ly[j-1]*lz[k])+(f[i-1][j-2][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j-2]*lz[k]))/((ly[j-2]+ly[j-1])*(lz[k-1]+lz[k]));
            double E_sb=((f[i+1][j-1][k]*ly[j-2]*lz[k-1])+(f[i+1][j-2][k-1]*ly[j-1]*lz[k])+(f[i+1][j-2][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j-2]*lz[k]))/((ly[j-2]+ly[j-1])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_et=((f[i+1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_wt=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_wt=((f[i][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_wt=((f[i][j-2][k+1]*lx[i-1]*lz[k])+(f[i-1][j-2][k]*lx[i]*lz[k+1])+(f[i-1][j-2][k+1]*lx[i]*lz[k])+(f[i][j-2][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_et=((f[i+1][j-2][k+1]*lx[i]*lz[k])+(f[i][j-2][k]*lx[i+1]*lz[k+1])+(f[i][j-2][k+1]*lx[i+1]*lz[k])+(f[i+1][j-2][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_eb=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_eb=((f[i+1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_wb=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_wb=((f[i][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_wb=((f[i][j-2][k]*lx[i-1]*lz[k-1])+(f[i-1][j-2][k-1]*lx[i]*lz[k])+(f[i-1][j-2][k]*lx[i]*lz[k-1])+(f[i][j-2][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_eb=((f[i+1][j-2][k]*lx[i]*lz[k-1])+(f[i][j-2][k-1]*lx[i+1]*lz[k])+(f[i][j-2][k]*lx[i+1]*lz[k-1])+(f[i+1][j-2][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j]+ly[j-1]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j-1]+ly[j-2]));
            double term_12=(P_et-S_et)/(0.5*(ly[j-1]+ly[j-2]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j]+ly[j-1]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j-1]+ly[j-2]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j-1]+ly[j-2]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_ne=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double T_nw=((f[i][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_nw=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double T_sw=((f[i][j-1][k+1]*lx[i-1]*ly[j-2])+(f[i-1][j-2][k+1]*lx[i]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j-2])+(f[i][j-2][k+1]*lx[i-1]*ly[j-1]))/((lx[i-1]+lx[i])*(ly[j-2]+ly[j-1]));
            double P_sw=((f[i][j-1][k]*lx[i-1]*ly[j-2])+(f[i-1][j-2][k]*lx[i]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j-2])+(f[i][j-2][k]*lx[i-1]*ly[j-1]))/((lx[i-1]+lx[i])*(ly[j-2]+ly[j-1]));
            double T_se=((f[i+1][j-1][k+1]*lx[i]*ly[j-2])+(f[i][j-2][k+1]*lx[i+1]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j-2])+(f[i+1][j-2][k+1]*lx[i]*ly[j-1]))/((lx[i]+lx[i+1])*(ly[j-2]+ly[j-1]));
            double P_se=((f[i+1][j-1][k]*lx[i]*ly[j-2])+(f[i][j-2][k]*lx[i+1]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j-2])+(f[i+1][j-2][k]*lx[i]*ly[j-1]))/((lx[i]+lx[i+1])*(ly[j-2]+ly[j-1]));
            double B_ne=((f[i+1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_nw=((f[i][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_sw=((f[i][j-1][k-1]*lx[i-1]*ly[j-2])+(f[i-1][j-2][k-1]*lx[i]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j-2])+(f[i][j-2][k-1]*lx[i-1]*ly[j-1]))/((lx[i-1]+lx[i])*(ly[j-2]+ly[j-1]));
            double B_se=((f[i+1][j-1][k-1]*lx[i]*ly[j-2])+(f[i][j-2][k-1]*lx[i+1]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j-2])+(f[i+1][j-2][k-1]*lx[i]*ly[j-1]))/((lx[i]+lx[i+1])*(ly[j-2]+ly[j-1]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(nory/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along +y

            stages_CICSAM();

            fs=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }
    }												//END Y-SWEEP...................................................

    else											//BEGIN Z-SWEEP.................................................
    {
        if (wt>0)										//fluxing is OUTWARD (wt is +ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j][k];
            f_ACCEPT=f[i][j][k+1];
            f_UPWIND=f[i][j][k-1];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j,k]
            double wt_D,wb_D;
            face_COUR=0.0;
            {
            wt_D=(lz[k+1]*w[i][j][k]+lz[k]*w[i][j][k+1])/(lz[k+1]+lz[k]);
            wb_D=(lz[k]*w[i][j][k-1]+lz[k-1]*w[i][j][k])/(lz[k]+lz[k-1]);
                face_COUR=fmax(wt_D*dt/lz[k],0.0)+fmax(-wb_D*dt/lz[k],0.0);
            }

            //calculation of the interface normal (Parker and Youngs' method)
            double E_nt=((f[i+1][j+1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_nt=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_nt=((f[i-1][j+1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_st=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_st=((f[i-1][j][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_st=((f[i+1][j][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_nb=((f[i+1][j+1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_nb=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_nb=((f[i-1][j+1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_sb=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_sb=((f[i-1][j][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_sb=((f[i+1][j][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_et=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_wt=((f[i][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_wt=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_wt=((f[i][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_et=((f[i+1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_eb=((f[i+1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_eb=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_wb=((f[i][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_wb=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_wb=((f[i][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_eb=((f[i+1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_ne=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double T_nw=((f[i][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_nw=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_sw=((f[i][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_sw=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double T_se=((f[i+1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_se=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_ne=((f[i+1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_nw=((f[i][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_sw=((f[i][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_se=((f[i+1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(norz/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along +z

            stages_CICSAM();

            ft=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }

        else											//fluxing is INWARD (wt is -ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j][k+1];
            f_ACCEPT=f[i][j][k];
            f_UPWIND=f[i][j][k+2];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j,k+1]
            double wt_D,wb_D;
            face_COUR=0.0;
            {
            wt_D=(lz[k+2]*w[i][j][k+1]+lz[k+1]*w[i][j][k+2])/(lz[k+2]+lz[k+1]);
            wb_D=(lz[k+1]*w[i][j][k]+lz[k]*w[i][j][k+1])/(lz[k+1]+lz[k]);
                face_COUR=fmax(wt_D*dt/lz[k],0.0)+fmax(-wb_D*dt/lz[k],0.0);
            }

            /*...calculation of the interface normal (Parker and Youngs' method) [stencil shifts one cell along +k]...*/
            double E_nt=((f[i+1][j+1][k+2]*ly[j]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k+2])+(f[i+1][j][k+2]*ly[j+1]*lz[k+1])+(f[i+1][j+1][k+1]*ly[j]*lz[k+2]))/((ly[j]+ly[j+1])*(lz[k+1]+lz[k+2]));
            double P_nt=((f[i][j+1][k+2]*ly[j]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k+2])+(f[i][j][k+2]*ly[j+1]*lz[k+1])+(f[i][j+1][k+1]*ly[j]*lz[k+2]))/((ly[j]+ly[j+1])*(lz[k+1]+lz[k+2]));
            double W_nt=((f[i-1][j+1][k+2]*ly[j]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k+2])+(f[i-1][j][k+2]*ly[j+1]*lz[k+1])+(f[i-1][j+1][k+1]*ly[j]*lz[k+2]))/((ly[j]+ly[j+1])*(lz[k+1]+lz[k+2]));
            double P_st=((f[i][j][k+2]*ly[j-1]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k+2])+(f[i][j-1][k+2]*ly[j]*lz[k+1])+(f[i][j][k+1]*ly[j-1]*lz[k+2]))/((ly[j-1]+ly[j])*(lz[k+1]+lz[k+2]));
            double W_st=((f[i-1][j][k+2]*ly[j-1]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k+2])+(f[i-1][j-1][k+2]*ly[j]*lz[k+1])+(f[i-1][j][k+1]*ly[j-1]*lz[k+2]))/((ly[j-1]+ly[j])*(lz[k+1]+lz[k+2]));
            double E_st=((f[i+1][j][k+2]*ly[j-1]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k+2])+(f[i+1][j-1][k+2]*ly[j]*lz[k+1])+(f[i+1][j][k+1]*ly[j-1]*lz[k+2]))/((ly[j-1]+ly[j])*(lz[k+1]+lz[k+2]));
            double E_nb=((f[i+1][j+1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_nb=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_nb=((f[i-1][j+1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_sb=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_sb=((f[i-1][j][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_sb=((f[i+1][j][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j+1][k+2]*lx[i]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k+2])+(f[i][j+1][k+2]*lx[i+1]*lz[k+1])+(f[i+1][j+1][k+1]*lx[i]*lz[k+2]))/((lx[i]+lx[i+1])*(lz[k+1]+lz[k+2]));
            double P_et=((f[i+1][j][k+2]*lx[i]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k+2])+(f[i][j][k+2]*lx[i+1]*lz[k+1])+(f[i+1][j][k+1]*lx[i]*lz[k+2]))/((lx[i]+lx[i+1])*(lz[k+1]+lz[k+2]));
            double N_wt=((f[i][j+1][k+2]*lx[i-1]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k+2])+(f[i-1][j+1][k+2]*lx[i]*lz[k+1])+(f[i][j+1][k+1]*lx[i-1]*lz[k+2]))/((lx[i-1]+lx[i])*(lz[k+1]+lz[k+2]));
            double P_wt=((f[i][j][k+2]*lx[i-1]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k+2])+(f[i-1][j][k+2]*lx[i]*lz[k+1])+(f[i][j][k+1]*lx[i-1]*lz[k+2]))/((lx[i-1]+lx[i])*(lz[k+1]+lz[k+2]));
            double S_wt=((f[i][j-1][k+2]*lx[i-1]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k+2])+(f[i-1][j-1][k+2]*lx[i]*lz[k+1])+(f[i][j-1][k+1]*lx[i-1]*lz[k+2]))/((lx[i-1]+lx[i])*(lz[k+1]+lz[k+2]));
            double S_et=((f[i+1][j-1][k+2]*lx[i]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k+2])+(f[i][j-1][k+2]*lx[i+1]*lz[k+1])+(f[i+1][j-1][k+1]*lx[i]*lz[k+2]))/((lx[i]+lx[i+1])*(lz[k+1]+lz[k+2]));
            double N_eb=((f[i+1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_eb=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_wb=((f[i][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_wb=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_wb=((f[i][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_eb=((f[i+1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j+1][k+2]*lx[i]*ly[j])+(f[i][j][k+2]*lx[i+1]*ly[j+1])+(f[i][j+1][k+2]*lx[i+1]*ly[j])+(f[i+1][j][k+2]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_ne=((f[i+1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double T_nw=((f[i][j+1][k+2]*lx[i-1]*ly[j])+(f[i-1][j][k+2]*lx[i]*ly[j+1])+(f[i-1][j+1][k+2]*lx[i]*ly[j])+(f[i][j][k+2]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_nw=((f[i][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_sw=((f[i][j][k+2]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+2]*lx[i]*ly[j])+(f[i-1][j][k+2]*lx[i]*ly[j-1])+(f[i][j-1][k+2]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_sw=((f[i][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double T_se=((f[i+1][j][k+2]*lx[i]*ly[j-1])+(f[i][j-1][k+2]*lx[i+1]*ly[j])+(f[i][j][k+2]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+2]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_se=((f[i+1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_ne=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_nw=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_sw=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_se=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(-norz/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along -z

            stages_CICSAM();

            ft=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }

        if (wb<0)										//fluxing is OUTWARD (wb is -ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j][k];
            f_ACCEPT=f[i][j][k-1];
            f_UPWIND=f[i][j][k+1];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j,k]
            double wt_D,wb_D;
            face_COUR=0.0;
            {
            wt_D=(lz[k+1]*w[i][j][k]+lz[k]*w[i][j][k+1])/(lz[k+1]+lz[k]);
            wb_D=(lz[k]*w[i][j][k-1]+lz[k-1]*w[i][j][k])/(lz[k]+lz[k-1]);
                face_COUR=fmax(wt_D*dt/lz[k],0.0)+fmax(-wb_D*dt/lz[k],0.0);
            }

            //calculation of the interface normal (Parker and Youngs' method)
            double E_nt=((f[i+1][j+1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k+1])+(f[i+1][j][k+1]*ly[j+1]*lz[k])+(f[i+1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_nt=((f[i][j+1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k+1])+(f[i][j][k+1]*ly[j+1]*lz[k])+(f[i][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double W_nt=((f[i-1][j+1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k+1])+(f[i-1][j][k+1]*ly[j+1]*lz[k])+(f[i-1][j+1][k]*ly[j]*lz[k+1]))/((ly[j]+ly[j+1])*(lz[k]+lz[k+1]));
            double P_st=((f[i][j][k+1]*ly[j-1]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k+1])+(f[i][j-1][k+1]*ly[j]*lz[k])+(f[i][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double W_st=((f[i-1][j][k+1]*ly[j-1]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k+1])+(f[i-1][j-1][k+1]*ly[j]*lz[k])+(f[i-1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_st=((f[i+1][j][k+1]*ly[j-1]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k+1])+(f[i+1][j-1][k+1]*ly[j]*lz[k])+(f[i+1][j][k]*ly[j-1]*lz[k+1]))/((ly[j-1]+ly[j])*(lz[k]+lz[k+1]));
            double E_nb=((f[i+1][j+1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_nb=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_nb=((f[i-1][j+1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_sb=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_sb=((f[i-1][j][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_sb=((f[i+1][j][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k+1])+(f[i][j+1][k+1]*lx[i+1]*lz[k])+(f[i+1][j+1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double P_et=((f[i+1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k+1])+(f[i][j][k+1]*lx[i+1]*lz[k])+(f[i+1][j][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_wt=((f[i][j+1][k+1]*lx[i-1]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k+1])+(f[i-1][j+1][k+1]*lx[i]*lz[k])+(f[i][j+1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double P_wt=((f[i][j][k+1]*lx[i-1]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k+1])+(f[i-1][j][k+1]*lx[i]*lz[k])+(f[i][j][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_wt=((f[i][j-1][k+1]*lx[i-1]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k+1])+(f[i-1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i-1]*lz[k+1]))/((lx[i-1]+lx[i])*(lz[k]+lz[k+1]));
            double S_et=((f[i+1][j-1][k+1]*lx[i]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k+1])+(f[i][j-1][k+1]*lx[i+1]*lz[k])+(f[i+1][j-1][k]*lx[i]*lz[k+1]))/((lx[i]+lx[i+1])*(lz[k]+lz[k+1]));
            double N_eb=((f[i+1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_eb=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_wb=((f[i][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_wb=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_wb=((f[i][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_eb=((f[i+1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j+1])+(f[i][j+1][k+1]*lx[i+1]*ly[j])+(f[i+1][j][k+1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_ne=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double T_nw=((f[i][j+1][k+1]*lx[i-1]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j+1])+(f[i-1][j+1][k+1]*lx[i]*ly[j])+(f[i][j][k+1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_nw=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_sw=((f[i][j][k+1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k+1]*lx[i]*ly[j])+(f[i-1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_sw=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double T_se=((f[i+1][j][k+1]*lx[i]*ly[j-1])+(f[i][j-1][k+1]*lx[i+1]*ly[j])+(f[i][j][k+1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k+1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_se=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_ne=((f[i+1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_nw=((f[i][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_sw=((f[i][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_se=((f[i+1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k+1]+lz[k]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k+1]+lz[k]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k+1]+lz[k]));
            double term_20=(T_se-P_se)/(0.5*(lz[k+1]+lz[k]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_24=(P_se-B_se)/(0.5*(lz[k]+lz[k-1]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(-norz/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along -z

            stages_CICSAM();

            fb=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }

        else											//fluxing is INWARD (wb is +ve)
        {
            //definition of DONOR, ACCEPTOR and UPWIND cells ; face Courant number ; interface angle
            f_DONOR=f[i][j][k-1];
            f_ACCEPT=f[i][j][k];
            f_UPWIND=f[i][j][k-2];
            ftild_DONOR=(f_DONOR-f_UPWIND)/(f_ACCEPT-f_UPWIND);
            //calculating face-Courant the right way  [DONOR is i,j,k-1]
            double wt_D,wb_D;
            face_COUR=0.0;
            {
            wt_D=(lz[k]*w[i][j][k-1]+lz[k-1]*w[i][j][k])/(lz[k]+lz[k-1]);
            wb_D=(lz[k-1]*w[i][j][k-2]+lz[k-2]*w[i][j][k-1])/(lz[k-1]+lz[k-2]);
                face_COUR=fmax(wt_D*dt/lz[k],0.0)+fmax(-wb_D*dt/lz[k],0.0);
            }

            /*...calculation of the interface normal (Parker and Youngs' method) [stencil shifts one cell along -k]...*/
            double E_nt=((f[i+1][j+1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k])+(f[i+1][j][k]*ly[j+1]*lz[k-1])+(f[i+1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_nt=((f[i][j+1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k])+(f[i][j][k]*ly[j+1]*lz[k-1])+(f[i][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double W_nt=((f[i-1][j+1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k])+(f[i-1][j][k]*ly[j+1]*lz[k-1])+(f[i-1][j+1][k-1]*ly[j]*lz[k]))/((ly[j]+ly[j+1])*(lz[k-1]+lz[k]));
            double P_st=((f[i][j][k]*ly[j-1]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k])+(f[i][j-1][k]*ly[j]*lz[k-1])+(f[i][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double W_st=((f[i-1][j][k]*ly[j-1]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k])+(f[i-1][j-1][k]*ly[j]*lz[k-1])+(f[i-1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_st=((f[i+1][j][k]*ly[j-1]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k])+(f[i+1][j-1][k]*ly[j]*lz[k-1])+(f[i+1][j][k-1]*ly[j-1]*lz[k]))/((ly[j-1]+ly[j])*(lz[k-1]+lz[k]));
            double E_nb=((f[i+1][j+1][k-1]*ly[j]*lz[k-2])+(f[i+1][j][k-2]*ly[j+1]*lz[k-1])+(f[i+1][j][k-1]*ly[j+1]*lz[k-2])+(f[i+1][j+1][k-2]*ly[j]*lz[k-1]))/((ly[j]+ly[j+1])*(lz[k-2]+lz[k-1]));
            double P_nb=((f[i][j+1][k-1]*ly[j]*lz[k-2])+(f[i][j][k-2]*ly[j+1]*lz[k-1])+(f[i][j][k-1]*ly[j+1]*lz[k-2])+(f[i][j+1][k-2]*ly[j]*lz[k-1]))/((ly[j]+ly[j+1])*(lz[k-2]+lz[k-1]));
            double W_nb=((f[i-1][j+1][k-1]*ly[j]*lz[k-2])+(f[i-1][j][k-2]*ly[j+1]*lz[k-1])+(f[i-1][j][k-1]*ly[j+1]*lz[k-2])+(f[i-1][j+1][k-2]*ly[j]*lz[k-1]))/((ly[j]+ly[j+1])*(lz[k-2]+lz[k-1]));
            double P_sb=((f[i][j][k-1]*ly[j-1]*lz[k-2])+(f[i][j-1][k-2]*ly[j]*lz[k-1])+(f[i][j-1][k-1]*ly[j]*lz[k-2])+(f[i][j][k-2]*ly[j-1]*lz[k-1]))/((ly[j-1]+ly[j])*(lz[k-2]+lz[k-1]));
            double W_sb=((f[i-1][j][k-1]*ly[j-1]*lz[k-2])+(f[i-1][j-1][k-2]*ly[j]*lz[k-1])+(f[i-1][j-1][k-1]*ly[j]*lz[k-2])+(f[i-1][j][k-2]*ly[j-1]*lz[k-1]))/((ly[j-1]+ly[j])*(lz[k-2]+lz[k-1]));
            double E_sb=((f[i+1][j][k-1]*ly[j-1]*lz[k-2])+(f[i+1][j-1][k-2]*ly[j]*lz[k-1])+(f[i+1][j-1][k-1]*ly[j]*lz[k-2])+(f[i+1][j][k-2]*ly[j-1]*lz[k-1]))/((ly[j-1]+ly[j])*(lz[k-2]+lz[k-1]));

            double term_1=(E_nt-P_nt)/(0.5*(lx[i+1]+lx[i]));
            double term_2=(P_nt-W_nt)/(0.5*(lx[i]+lx[i-1]));
            double term_3=(P_st-W_st)/(0.5*(lx[i]+lx[i-1]));
            double term_4=(E_st-P_st)/(0.5*(lx[i+1]+lx[i]));
            double term_5=(E_nb-P_nb)/(0.5*(lx[i+1]+lx[i]));
            double term_6=(P_nb-W_nb)/(0.5*(lx[i]+lx[i-1]));
            double term_7=(P_sb-W_sb)/(0.5*(lx[i]+lx[i-1]));
            double term_8=(E_sb-P_sb)/(0.5*(lx[i+1]+lx[i]));

            norx=0.5*(term_1+term_2+term_3+term_4+term_5+term_6+term_7+term_8); 		/*X-COMPONENT OF THE INTERFACE NORMAL*/

            double N_et=((f[i+1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k])+(f[i][j+1][k]*lx[i+1]*lz[k-1])+(f[i+1][j+1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double P_et=((f[i+1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k])+(f[i][j][k]*lx[i+1]*lz[k-1])+(f[i+1][j][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_wt=((f[i][j+1][k]*lx[i-1]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k])+(f[i-1][j+1][k]*lx[i]*lz[k-1])+(f[i][j+1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double P_wt=((f[i][j][k]*lx[i-1]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k])+(f[i-1][j][k]*lx[i]*lz[k-1])+(f[i][j][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_wt=((f[i][j-1][k]*lx[i-1]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k])+(f[i-1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i-1]*lz[k]))/((lx[i-1]+lx[i])*(lz[k-1]+lz[k]));
            double S_et=((f[i+1][j-1][k]*lx[i]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k])+(f[i][j-1][k]*lx[i+1]*lz[k-1])+(f[i+1][j-1][k-1]*lx[i]*lz[k]))/((lx[i]+lx[i+1])*(lz[k-1]+lz[k]));
            double N_eb=((f[i+1][j+1][k-1]*lx[i]*lz[k-2])+(f[i][j+1][k-2]*lx[i+1]*lz[k-1])+(f[i][j+1][k-1]*lx[i+1]*lz[k-2])+(f[i+1][j+1][k-2]*lx[i]*lz[k-1]))/((lx[i]+lx[i+1])*(lz[k-2]+lz[k-1]));
            double P_eb=((f[i+1][j][k-1]*lx[i]*lz[k-2])+(f[i][j][k-2]*lx[i+1]*lz[k-1])+(f[i][j][k-1]*lx[i+1]*lz[k-2])+(f[i+1][j][k-2]*lx[i]*lz[k-1]))/((lx[i]+lx[i+1])*(lz[k-2]+lz[k-1]));
            double N_wb=((f[i][j+1][k-1]*lx[i-1]*lz[k-2])+(f[i-1][j+1][k-2]*lx[i]*lz[k-1])+(f[i-1][j+1][k-1]*lx[i]*lz[k-2])+(f[i][j+1][k-2]*lx[i-1]*lz[k-1]))/((lx[i-1]+lx[i])*(lz[k-2]+lz[k-1]));
            double P_wb=((f[i][j][k-1]*lx[i-1]*lz[k-2])+(f[i-1][j][k-2]*lx[i]*lz[k-1])+(f[i-1][j][k-1]*lx[i]*lz[k-2])+(f[i][j][k-2]*lx[i-1]*lz[k-1]))/((lx[i-1]+lx[i])*(lz[k-2]+lz[k-1]));
            double S_wb=((f[i][j-1][k-1]*lx[i-1]*lz[k-2])+(f[i-1][j-1][k-2]*lx[i]*lz[k-1])+(f[i-1][j-1][k-1]*lx[i]*lz[k-2])+(f[i][j-1][k-2]*lx[i-1]*lz[k-1]))/((lx[i-1]+lx[i])*(lz[k-2]+lz[k-1]));
            double S_eb=((f[i+1][j-1][k-1]*lx[i]*lz[k-2])+(f[i][j-1][k-2]*lx[i+1]*lz[k-1])+(f[i][j-1][k-1]*lx[i+1]*lz[k-2])+(f[i+1][j-1][k-2]*lx[i]*lz[k-1]))/((lx[i]+lx[i+1])*(lz[k-2]+lz[k-1]));

            double term_9=(N_et-P_et)/(0.5*(ly[j+1]+ly[j]));
            double term_10=(N_wt-P_wt)/(0.5*(ly[j+1]+ly[j]));
            double term_11=(P_wt-S_wt)/(0.5*(ly[j]+ly[j-1]));
            double term_12=(P_et-S_et)/(0.5*(ly[j]+ly[j-1]));
            double term_13=(N_eb-P_eb)/(0.5*(ly[j+1]+ly[j]));
            double term_14=(N_wb-P_wb)/(0.5*(ly[j+1]+ly[j]));
            double term_15=(P_wb-S_wb)/(0.5*(ly[j]+ly[j-1]));
            double term_16=(P_eb-S_eb)/(0.5*(ly[j]+ly[j-1]));

            nory=0.5*(term_9+term_10+term_11+term_12+term_13+term_14+term_15+term_16);  	/*Y-COMPONENT OF THE INTERFACE NORMAL*/

            double T_ne=((f[i+1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j+1])+(f[i][j+1][k]*lx[i+1]*ly[j])+(f[i+1][j][k]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double P_ne=((f[i+1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j+1])+(f[i][j+1][k-1]*lx[i+1]*ly[j])+(f[i+1][j][k-1]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double T_nw=((f[i][j+1][k]*lx[i-1]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j+1])+(f[i-1][j+1][k]*lx[i]*ly[j])+(f[i][j][k]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double P_nw=((f[i][j+1][k-1]*lx[i-1]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j+1])+(f[i-1][j+1][k-1]*lx[i]*ly[j])+(f[i][j][k-1]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double T_sw=((f[i][j][k]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k]*lx[i]*ly[j])+(f[i-1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double P_sw=((f[i][j][k-1]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-1]*lx[i]*ly[j])+(f[i-1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double T_se=((f[i+1][j][k]*lx[i]*ly[j-1])+(f[i][j-1][k]*lx[i+1]*ly[j])+(f[i][j][k]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double P_se=((f[i+1][j][k-1]*lx[i]*ly[j-1])+(f[i][j-1][k-1]*lx[i+1]*ly[j])+(f[i][j][k-1]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-1]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));
            double B_ne=((f[i+1][j+1][k-2]*lx[i]*ly[j])+(f[i][j][k-2]*lx[i+1]*ly[j+1])+(f[i][j+1][k-2]*lx[i+1]*ly[j])+(f[i+1][j][k-2]*lx[i]*ly[j+1]))/((lx[i]+lx[i+1])*(ly[j]+ly[j+1]));
            double B_nw=((f[i][j+1][k-2]*lx[i-1]*ly[j])+(f[i-1][j][k-2]*lx[i]*ly[j+1])+(f[i-1][j+1][k-2]*lx[i]*ly[j])+(f[i][j][k-2]*lx[i-1]*ly[j+1]))/((lx[i-1]+lx[i])*(ly[j]+ly[j+1]));
            double B_sw=((f[i][j][k-2]*lx[i-1]*ly[j-1])+(f[i-1][j-1][k-2]*lx[i]*ly[j])+(f[i-1][j][k-2]*lx[i]*ly[j-1])+(f[i][j-1][k-2]*lx[i-1]*ly[j]))/((lx[i-1]+lx[i])*(ly[j-1]+ly[j]));
            double B_se=((f[i+1][j][k-2]*lx[i]*ly[j-1])+(f[i][j-1][k-2]*lx[i+1]*ly[j])+(f[i][j][k-2]*lx[i+1]*ly[j-1])+(f[i+1][j-1][k-2]*lx[i]*ly[j]))/((lx[i]+lx[i+1])*(ly[j-1]+ly[j]));

            double term_17=(T_ne-P_ne)/(0.5*(lz[k]+lz[k-1]));
            double term_18=(T_nw-P_nw)/(0.5*(lz[k]+lz[k-1]));
            double term_19=(T_sw-P_sw)/(0.5*(lz[k]+lz[k-1]));
            double term_20=(T_se-P_se)/(0.5*(lz[k]+lz[k-1]));
            double term_21=(P_ne-B_ne)/(0.5*(lz[k-1]+lz[k-2]));
            double term_22=(P_nw-B_nw)/(0.5*(lz[k-1]+lz[k-2]));
            double term_23=(P_sw-B_sw)/(0.5*(lz[k-1]+lz[k-2]));
            double term_24=(P_se-B_se)/(0.5*(lz[k-1]+lz[k-2]));

            norz=0.5*(term_17+term_18+term_19+term_20+term_21+term_22+term_23+term_24); 	/*Z-COMPONENT OF THE INTERFACE NORMAL*/

            theta_f=acos(norz/sqrt(norx*norx+nory*nory+norz*norz));				//unit vector connecting D and A is along +z

            stages_CICSAM();

            fb=(1.0-beta_f)*f_DONOR+beta_f*f_ACCEPT;
        }
    }												//END Z-SWEEP...................................................
}

/********************************/								/**************************************************************/
/*...REDISTRIBUTION ALGORITHM...*/								/*...Saincher & Banerjee (2015) - Numerical Heat Transfer-B...*/
/********************************/								/**************************************************************/
void retribution()
{
do
{
forlorn=0;
for (i=IMG+1;i<=NX+IMG;i++)
{
    for (j=IMG+1;j<=NY+IMG;j++)
    {
        for (k=IMG+1;k<=NZ+IMG;k++)
        {
            if (f[i][j][k]<0.0-epsilon)								//undershoot
            {
            usc=usc+1;
            int index=1;
             do{
                min=10.0;
                for (int l=i-index;l<=i+index;l++)
                {
                    for (int m=j-index;m<=j+index;m++)
                    {
                        for (int n=k-index;n<=k+index;n++)
                        {
                            if ((l==i&&m==j&&n==k)||l<(IMG+1)||m<(IMG+1)||n<(IMG+1)||l>(NX+IMG)||m>(NY+IMG)||n>(NZ+IMG))//skipping imaginary cells
                            continue;
                            else
                            {
                            if (f[l][m][n]<=min&&f[l][m][n]!=0.0)
                            {
                            min=f[l][m][n];
                            a1=l;
                            b1=m;
                            c1=n;
                            }
                            }
                        }
                    }
                }
                if (f[a1][b1][c1]==0.0)
                {
                index=index+1;
                }
              temp=f[a1][b1][c1]+f[i][j][k];
                {
                if (temp>0.0)
                f[a1][b1][c1]=temp;
                else
                f[a1][b1][c1]=0.0;
                }
                {
                if (temp<0.0)
                f[i][j][k]=temp;
                else
                f[i][j][k]=0.0;
                }
               }while(f[i][j][k]<0.0-epsilon);
            }											//end UNDERSHOOT condition
            else if (f[i][j][k]>1.0+epsilon)							//overshoot
            {
            osc=osc+1;
            int index=1;
             do{
                max=-10.0;
                for (int l=i-index;l<=i+index;l++)
                {
                    for (int m=j-index;m<=j+index;m++)
                    {
                        for (int n=k-index;n<=k+index;n++)
                        {
                            if ((l==i&&m==j&&n==k)||l<(IMG+1)||m<(IMG+1)||n<(IMG+1)||l>(NX+IMG)||m>(NY+IMG)||n>(NZ+IMG))//skipping imaginary cells
                            continue;
                            else
                            {
                            if (f[l][m][n]>=max&&f[l][m][n]!=1.0)
                            {
                            max=f[l][m][n];
                            a1=l;
                            b1=m;
                            c1=n;
                            }
                            }
                        }
                    }
                }
                if (f[a1][b1][c1]==1.0)
                {
                index=index+1;
                }
              temp=f[a1][b1][c1]+f[i][j][k]-1.0;
                {
                if (temp<1.0)
                f[a1][b1][c1]=temp;
                else
                f[a1][b1][c1]=1.0;
                }
                {
                if (temp>1.0)
                f[i][j][k]=temp;
                else
                f[i][j][k]=1.0;
                }
               }while(f[i][j][k]>1.0+epsilon);
            }											//end OVERSHOOT condition
            else
            continue;
        }											//end k-LOOP
    }												//end j-LOOP
}												//end i-LOOP
for (i=IMG+1;i<=NX+IMG;i++)
{
    for (j=IMG+1;j<=NY+IMG;j++)
    {
        for (k=IMG+1;k<=NZ+IMG;k++)
        {
            if (f[i][j][k]>1.0+epsilon||f[i][j][k]<0.0-epsilon)
            {
            forlorn=1;
            printf("%d\t%d\t%d\t%f\n",i,j,k,f[i][j][k]);
            break;
            }
            else
            continue;
        }
    }
}
}while (forlorn>0);
}

//x-sweep\\...
void xsweep()
{
swp=0;
    for (i=IMG+1;i<=NX+IMG;i++)
    {
    for (j=IMG+1;j<=NY+IMG;j++)
    {
    for (k=IMG+1;k<=NZ+IMG;k++)
    {
        //estimation X-momentum, advecting VOF (collocated), from the underlying Staggered mesh
        ue=u[i+1][j][k];
        uw=u[i][j][k];

        f_CICSAM();

        fnew[i][j][k]=f[i][j][k]+(dt/lx[i])*(uw*fw-ue*fe)+c[i][j][k]*(dt/lx[i])*(ue-uw);    	/*...WEYMOUTH and YUE (2010)...*/
    }												//end k-loop
    }												//end j-loop
    }												//end i-loop
}

//y-sweep\\...
void ysweep()
{
swp=1;
    for (i=IMG+1;i<=NX+IMG;i++)
    {
    for (j=IMG+1;j<=NY+IMG;j++)
    {
    for (k=IMG+1;k<=NZ+IMG;k++)
    {
        //estimation Y-momentum, advecting VOF (collocated), from the underlying Staggered mesh
        vn=v[i][j+1][k];
        vs=v[i][j][k];

        f_CICSAM();

        fnew[i][j][k]=f[i][j][k]+(dt/ly[j])*(vs*fs-vn*fn)+c[i][j][k]*(dt/ly[j])*(vn-vs);    	/*...WEYMOUTH and YUE (2010)...*/
    }												//end k-loop
    }												//end j-loop
    }												//end i-loop
}

//z-sweep\\...
void zsweep()
{
swp=2;
    for (i=IMG+1;i<=NX+IMG;i++)
    {
    for (j=IMG+1;j<=NY+IMG;j++)
    {
    for (k=IMG+1;k<=NZ+IMG;k++)
    {
        //estimation Z-momentum, advecting VOF (collocated), from the underlying Staggered mesh
        wt=w[i][j][k+1];
        wb=w[i][j][k];

        f_CICSAM();

        fnew[i][j][k]=f[i][j][k]+(dt/lz[k])*(wb*fb-wt*ft)+c[i][j][k]*(dt/lz[k])*(wt-wb);    	/*...WEYMOUTH and YUE (2010)...*/
    }												//end k-loop
    }												//end j-loop
    }												//end i-loop
}

void vof_c_update()										//the "old" and "new" general volume fractions are the same after vof_c_update
{
    for (i=IMG+1;i<=NX+IMG;i++)
        for (j=IMG+1;j<=NY+IMG;j++)
            for (k=IMG+1;k<=NZ+IMG;k++)
                f[i][j][k]=fnew[i][j][k];
}

/*...FUNCTIONS NECESSARY FOR U-ADVECTION ROUTINE...*/
void uFOU()											//FIRST ORDER UPWIND (FOU)
{
    uadvel=phiuU;
}

void uSOU()											//SECOND ORDER UPWIND (SOU)
{
    uadvel=((2.0*VOLU+VOLUU)/(VOLU+VOLUU))*phiuU-(VOLU/(VOLU+VOLUU))*phiuUU;
}

void uQUICK()											//QUADRATIC UPSTREAM INTERPOLATION OF CONVECTION KINETICS (QUICK)
{
    uadvel=((VOLU*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLD+2.0*VOLU+VOLUU)))*phiuD+((VOLD*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLU+VOLUU)))*phiuU-((VOLD*VOLU)/((VOLU+VOLUU)*(VOLD+2.0*VOLU+VOLUU)))*phiuUU;
}

void uTVD()											//TOTAL VARIATION DIMINISHING (TVD)
{
    double xi=0.0;
    {
    xi=fabs((phiuU-phiuUU)/(phiuD-phiuUU));
        if (xi==0.0||xi==1.0)
        uadvel=phiuU;
        else if (xi>0.0&&xi<=0.3)
        uadvel=((3.0*VOLU+VOLUU)/(VOLU+VOLUU))*phiuU-((2.0*VOLU)/(VOLU+VOLUU))*phiuUU;
        else if (xi>0.3&&xi<=0.833)
        uadvel=((VOLU*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLD+2.0*VOLU+VOLUU)))*phiuD+((VOLD*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLU+VOLUU)))*phiuU-((VOLD*VOLU)/((VOLU+VOLUU)*(VOLD+2.0*VOLU+VOLUU)))*phiuUU;
        else if (xi>0.833&&xi<1.0)
        uadvel=phiuD;
        else
        uadvel=phiuU;
    }
}

void uFiOU()											//FIFTH ORDER UPWIND (FiOU)
{
    //cell-center distances from the face
    double xU1=-0.5*VOLU;									//upwind cell-center
    double xU2=-(VOLU+0.5*VOLUU);								//up-upwind cell-center
    double xU3=-(VOLU+VOLUU+0.5*VOLUUU);							//up-up-upwind cell-center
    double xD1=0.5*VOLD;									//downwind cell-center
    double xD2=VOLD+0.5*VOLDD;									//down-downwind cell-center

    uadvel=(1.0-((xU1*xU3*xD1*xD2)/((xU2-xU1)*(xU2-xU3)*(xU2-xD1)*(xU2-xD2)))-((xU1*xU2*xD1*xD2)/((xU3-xU1)*(xU3-xU2)*(xU3-xD1)*(xU3-xD2)))-((xU1*xU2*xU3*xD2)/((xD1-xU1)*(xD1-xU2)*(xD1-xU3)*(xD1-xD2)))-((xU1*xU2*xU3*xD1)/((xD2-xU1)*(xD2-xU2)*(xD2-xU3)*(xD2-xD1))))*phiuU+((xU1*xU3*xD1*xD2)/((xU2-xU1)*(xU2-xU3)*(xU2-xD1)*(xU2-xD2)))*phiuUU+((xU1*xU2*xD1*xD2)/((xU3-xU1)*(xU3-xU2)*(xU3-xD1)*(xU3-xD2)))*phiuUUU+((xU1*xU2*xU3*xD2)/((xD1-xU1)*(xD1-xU2)*(xD1-xU3)*(xD1-xD2)))*phiuD+((xU1*xU2*xU3*xD1)/((xD2-xU1)*(xD2-xU2)*(xD2-xU3)*(xD2-xD1)))*phiuDD;
}

void uSCHEMES()											//sch=4 is CD2
{
    if (sch==1)
    uFOU();
    else if (sch==2)
    uSOU();
    else if (sch==3)
    uQUICK();
    else if (sch==5)
    uTVD();
    else
    uFiOU();
}

/*...ADVECTION ROUTINE FOR U-VELOCITY...*/
void uadrou()
{
    /*...X-advection of U...*/
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            uadvel=ue;
        }
    else
        {
        if (ue>0)
        {
            phiuU=u[i][j][k];       	VOLU=0.5*(lx[i]+lx[i-1]);   				//staggered
            phiuUU=u[i-1][j][k];        VOLUU=0.5*(lx[i-1]+lx[i-2]);    			//staggered
            phiuUUU=u[i-2][j][k];       VOLUUU=0.5*(lx[i-2]+lx[i-3]);   			//staggered
            phiuD=u[i+1][j][k];         VOLD=0.5*(lx[i+1]+lx[i]);   				//staggered
            phiuDD=u[i+2][j][k];        VOLDD=0.5*(lx[i+2]+lx[i+1]);    			//staggered
            uSCHEMES();
        }
        else
        {
            phiuU=u[i+1][j][k];         VOLU=0.5*(lx[i+1]+lx[i]);   				//staggered 
            phiuUU=u[i+2][j][k];        VOLUU=0.5*(lx[i+2]+lx[i+1]);    			//staggered
            phiuUUU=u[i+3][j][k];       VOLUUU=0.5*(lx[i+3]+lx[i+2]);   			//staggered
            phiuD=u[i][j][k];       	VOLD=0.5*(lx[i]+lx[i-1]);   				//staggered
            phiuDD=u[i-1][j][k];        VOLDD=0.5*(lx[i-1]+lx[i-2]);    			//staggered
            uSCHEMES();
        }
        }
    uea=uadvel;
    }
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            uadvel=uw;
        }
    else
        {
        if (uw>0)
        {
            phiuU=u[i-1][j][k];         VOLU=0.5*(lx[i-1]+lx[i-2]); 				//staggered
            phiuUU=u[i-2][j][k];        VOLUU=0.5*(lx[i-2]+lx[i-3]);    			//staggered
            phiuUUU=u[i-3][j][k];       VOLUUU=0.5*(lx[i-3]+lx[i-4]);   			//staggered
            phiuD=u[i][j][k];       	VOLD=0.5*(lx[i]+lx[i-1]);   				//staggered
            phiuDD=u[i+1][j][k];        VOLDD=0.5*(lx[i+1]+lx[i]);  				//staggered
            uSCHEMES();
        }
        else
        {
            phiuU=u[i][j][k];       	VOLU=0.5*(lx[i]+lx[i-1]);   				//staggered
            phiuUU=u[i+1][j][k];        VOLUU=0.5*(lx[i+1]+lx[i]);  				//staggered
            phiuUUU=u[i+2][j][k];       VOLUUU=0.5*(lx[i+2]+lx[i+1]);   			//staggered
            phiuD=u[i-1][j][k];         VOLD=0.5*(lx[i-1]+lx[i-2]); 				//staggered
            phiuDD=u[i-2][j][k];        VOLDD=0.5*(lx[i-2]+lx[i-3]);    			//staggered
            uSCHEMES();
        }
        }
    uwa=uadvel;
    }
    /*...Y-advection of U...*/
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            uadvel=un;
        }
    else
        {
        if (vn>0)
        {
            phiuU=u[i][j][k];       	VOLU=ly[j]; 						//collocated
            phiuUU=u[i][j-1][k];        VOLUU=ly[j-1];  					//collocated
            phiuUUU=u[i][j-2][k];       VOLUUU=ly[j-2]; 					//collocated
            phiuD=u[i][j+1][k];         VOLD=ly[j+1];   					//collocated
            phiuDD=u[i][j+2][k];        VOLDD=ly[j+2];  					//collocated
            uSCHEMES();
        }
        else
        {
            phiuU=u[i][j+1][k];         VOLU=ly[j+1];   					//collocated
            phiuUU=u[i][j+2][k];        VOLUU=ly[j+2];  					//collocated
            phiuUUU=u[i][j+3][k];       VOLUUU=ly[j+3]; 					//collocated
            phiuD=u[i][j][k];       	VOLD=ly[j]; 						//collocated
            phiuDD=u[i][j-1][k];        VOLDD=ly[j-1];  					//collocated
            uSCHEMES();
        }
        }
    una=uadvel;
    }
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            uadvel=us;
        }
    else
        {
        if (vs>0)
        {
            phiuU=u[i][j-1][k];         VOLU=ly[j-1];   					//collocated
            phiuUU=u[i][j-2][k];        VOLUU=ly[j-2];  					//collocated
            phiuUUU=u[i][j-3][k];       VOLUUU=ly[j-3]; 					//collocated
            phiuD=u[i][j][k];       	VOLD=ly[j]; 						//collocated
            phiuDD=u[i][j+1][k];        VOLDD=ly[j+1];  					//collocated
            uSCHEMES();
        }
        else
        {
            phiuU=u[i][j][k];       	VOLU=ly[j]; 						//collocated
            phiuUU=u[i][j+1][k];        VOLUU=ly[j+1];  					//collocated
            phiuUUU=u[i][j+2][k];       VOLUUU=ly[j+2]; 					//collocated
            phiuD=u[i][j-1][k];         VOLD=ly[j-1];   					//collocated
            phiuDD=u[i][j-2][k];        VOLDD=ly[j-2];  					//collocated
            uSCHEMES();
        }
        }
    usa=uadvel;
    }
    /*...Z-advection of U...*/
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            uadvel=ut;
        }
    else
        {
        if (wt>0)
        {
            phiuU=u[i][j][k];       	VOLU=lz[k]; 						//collocated
            phiuUU=u[i][j][k-1];        VOLUU=lz[k-1];  					//collocated
            phiuUUU=u[i][j][k-2];       VOLUUU=lz[k-2]; 					//collocated
            phiuD=u[i][j][k+1];         VOLD=lz[k+1];   					//collocated
            phiuDD=u[i][j][k+2];        VOLDD=lz[k+2];  					//collocated
            uSCHEMES();
        }
        else
        {
            phiuU=u[i][j][k+1];         VOLU=lz[k+1];   					//collocated
            phiuUU=u[i][j][k+2];        VOLUU=lz[k+2];  					//collocated
            phiuUUU=u[i][j][k+3];       VOLUUU=lz[k+3]; 					//collocated
            phiuD=u[i][j][k];       	VOLD=lz[k]; 						//collocated
            phiuDD=u[i][j][k-1];        VOLDD=lz[k-1];  					//collocated
            uSCHEMES();
        }
        }
    uta=uadvel;
    }
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            uadvel=ub;
        }
    else
        {
        if (wb>0)
        {
            phiuU=u[i][j][k-1];         VOLU=lz[k-1];   					//collocated
            phiuUU=u[i][j][k-2];        VOLUU=lz[k-2];  					//collocated
            phiuUUU=u[i][j][k-3];       VOLUUU=lz[k-3]; 					//collocated
            phiuD=u[i][j][k];       	VOLD=lz[k]; 						//collocated
            phiuDD=u[i][j][k+1];        VOLDD=lz[k+1];  					//collocated
            uSCHEMES();
        }
        else
        {
            phiuU=u[i][j][k];       	VOLU=lz[k]; 						//collocated
            phiuUU=u[i][j][k+1];        VOLUU=lz[k+1];  					//collocated
            phiuUUU=u[i][j][k+2];       VOLUUU=lz[k+2]; 					//collocated
            phiuD=u[i][j][k-1];         VOLD=lz[k-1];   					//collocated
            phiuDD=u[i][j][k-2];        VOLDD=lz[k-2];  					//collocated
            uSCHEMES();
        }
        }
    uba=uadvel;
    }
}

/*...FUNCTIONS NECESSARY FOR V-ADVECTION ROUTINE...*/
void vFOU()											//FIRST ORDER UPWIND (FOU)
{
    vadvel=phivU;
}

void vSOU()											//SECOND ORDER UPWIND (SOU)
{
    vadvel=((2.0*VOLU+VOLUU)/(VOLU+VOLUU))*phivU-(VOLU/(VOLU+VOLUU))*phivUU;
}

void vQUICK()											//QUADRATIC UPSTREAM INTERPOLATION OF CONVECTION KINETICS (QUICK)
{
    vadvel=((VOLU*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLD+2.0*VOLU+VOLUU)))*phivD+((VOLD*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLU+VOLUU)))*phivU-((VOLD*VOLU)/((VOLU+VOLUU)*(VOLD+2.0*VOLU+VOLUU)))*phivUU;
}

void vTVD()											//TOTAL VARIATION DIMINISHING (TVD)
{
    double xi=0.0;
    {
    xi=fabs((phivU-phivUU)/(phivD-phivUU));
        if (xi==0.0||xi==1.0)
        vadvel=phivU;
        else if (xi>0.0&&xi<=0.3)
        vadvel=((3.0*VOLU+VOLUU)/(VOLU+VOLUU))*phivU-((2.0*VOLU)/(VOLU+VOLUU))*phivUU;
        else if (xi>0.3&&xi<=0.833)
        vadvel=((VOLU*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLD+2.0*VOLU+VOLUU)))*phivD+((VOLD*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLU+VOLUU)))*phivU-((VOLD*VOLU)/((VOLU+VOLUU)*(VOLD+2.0*VOLU+VOLUU)))*phivUU;
        else if (xi>0.833&&xi<1.0)
        vadvel=phivD;
        else
        vadvel=phivU;
    }
}

void vFiOU()											//FIFTH ORDER UPWIND (FiOU)
{
    //cell-center distances from the face
    double xU1=-0.5*VOLU;									//upwind cell-center
    double xU2=-(VOLU+0.5*VOLUU);								//up-upwind cell-center
    double xU3=-(VOLU+VOLUU+0.5*VOLUUU);							//up-up-upwind cell-center
    double xD1=0.5*VOLD;									//downwind cell-center
    double xD2=VOLD+0.5*VOLDD;									//down-downwind cell-center

    vadvel=(1.0-((xU1*xU3*xD1*xD2)/((xU2-xU1)*(xU2-xU3)*(xU2-xD1)*(xU2-xD2)))-((xU1*xU2*xD1*xD2)/((xU3-xU1)*(xU3-xU2)*(xU3-xD1)*(xU3-xD2)))-((xU1*xU2*xU3*xD2)/((xD1-xU1)*(xD1-xU2)*(xD1-xU3)*(xD1-xD2)))-((xU1*xU2*xU3*xD1)/((xD2-xU1)*(xD2-xU2)*(xD2-xU3)*(xD2-xD1))))*phivU+((xU1*xU3*xD1*xD2)/((xU2-xU1)*(xU2-xU3)*(xU2-xD1)*(xU2-xD2)))*phivUU+((xU1*xU2*xD1*xD2)/((xU3-xU1)*(xU3-xU2)*(xU3-xD1)*(xU3-xD2)))*phivUUU+((xU1*xU2*xU3*xD2)/((xD1-xU1)*(xD1-xU2)*(xD1-xU3)*(xD1-xD2)))*phivD+((xU1*xU2*xU3*xD1)/((xD2-xU1)*(xD2-xU2)*(xD2-xU3)*(xD2-xD1)))*phivDD;
}

void vSCHEMES()											//sch=4 is CD2
{
    if (sch==1)
    vFOU();
    else if (sch==2)
    vSOU();
    else if (sch==3)
    vQUICK();
    else if (sch==5)
    vTVD();
    else
    vFiOU();
}

/*...ADVECTION ROUTINE FOR V-VELOCITY...*/
void vadrou()
{
    /*...X-advection of V...*/
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            vadvel=ve;
        }
    else
        {
        if (ue>0)
        {
            phivU=v[i][j][k];       	VOLU=lx[i]; 						//collocated
            phivUU=v[i-1][j][k];        VOLUU=lx[i-1];  					//collocated
            phivUUU=v[i-2][j][k];       VOLUUU=lx[i-2]; 					//collocated
            phivD=v[i+1][j][k];         VOLD=lx[i+1];   					//collocated
            phivDD=v[i+2][j][k];        VOLDD=lx[i+2];  					//collocated
            vSCHEMES();
        }
        else
        {
            phivU=v[i+1][j][k];         VOLU=lx[i+1];   					//collocated
            phivUU=v[i+2][j][k];        VOLUU=lx[i+2];  					//collocated
            phivUUU=v[i+3][j][k];       VOLUUU=lx[i+3]; 					//collocated
            phivD=v[i][j][k];       	VOLD=lx[i]; 						//collocated
            phivDD=v[i-1][j][k];        VOLDD=lx[i-1];  					//collocated
            vSCHEMES();
        }
        }
    vea=vadvel;
    }
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            vadvel=vw;
        }
    else
        {
        if (uw>0)
        {
            phivU=v[i-1][j][k];         VOLU=lx[i-1];   					//collocated
            phivUU=v[i-2][j][k];        VOLUU=lx[i-2];  					//collocated
            phivUUU=v[i-3][j][k];       VOLUUU=lx[i-3]; 					//collocated
            phivD=v[i][j][k];       	VOLD=lx[i]; 						//collocated
            phivDD=v[i+1][j][k];        VOLDD=lx[i+1];  					//collocated
            vSCHEMES();
        }
        else
        {
            phivU=v[i][j][k];       	VOLU=lx[i]; 						//collocated
            phivUU=v[i+1][j][k];        VOLUU=lx[i+1];  					//collocated
            phivUUU=v[i+2][j][k];       VOLUUU=lx[i+2]; 					//collocated
            phivD=v[i-1][j][k];         VOLD=lx[i-1];   					//collocated
            phivDD=v[i-2][j][k];        VOLDD=lx[i-2];  					//collocated
            vSCHEMES();
        }
        }
    vwa=vadvel;
    }
    /*...Y-advection of V...*/
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            vadvel=vn;
        }
    else
        {
        if (vn>0)
        {
            phivU=v[i][j][k];       	VOLU=0.5*(ly[j]+ly[j-1]);   				//staggered
            phivUU=v[i][j-1][k];        VOLUU=0.5*(ly[j-1]+ly[j-2]);    			//staggered
            phivUUU=v[i][j-2][k];       VOLUUU=0.5*(ly[j-2]+ly[j-3]);   			//staggered
            phivD=v[i][j+1][k];         VOLD=0.5*(ly[j+1]+ly[j]);   				//staggered
            phivDD=v[i][j+2][k];        VOLDD=0.5*(ly[j+2]+ly[j+1]);    			//staggered
            vSCHEMES();
        }
        else
        {
            phivU=v[i][j+1][k];         VOLU=0.5*(ly[j+1]+ly[j]);   				//staggered
            phivUU=v[i][j+2][k];        VOLUU=0.5*(ly[j+2]+ly[j+1]);    			//staggered
            phivUUU=v[i][j+3][k];       VOLUUU=0.5*(ly[j+3]+ly[j+2]);   			//staggered
            phivD=v[i][j][k];       	VOLD=0.5*(ly[j]+ly[j-1]);   				//staggered
            phivDD=v[i][j-1][k];        VOLDD=0.5*(ly[j-1]+ly[j-2]);    			//staggered
            vSCHEMES();
        }
        }
    vna=vadvel;
    }
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            vadvel=vs;
        }
    else
        {
        if (vs>0)
        {
            phivU=v[i][j-1][k];         VOLU=0.5*(ly[j-1]+ly[j-2]); 				//staggered
            phivUU=v[i][j-2][k];        VOLUU=0.5*(ly[j-2]+ly[j-3]);    			//staggered
            phivUUU=v[i][j-3][k];       VOLUUU=0.5*(ly[j-3]+ly[j-4]);   			//staggered
            phivD=v[i][j][k];       	VOLD=0.5*(ly[j]+ly[j-1]);   				//staggered
            phivDD=v[i][j+1][k];        VOLDD=0.5*(ly[j+1]+ly[j]);  				//staggered
            vSCHEMES();
        }
        else
        {
            phivU=v[i][j][k];       	VOLU=0.5*(ly[j]+ly[j-1]);   				//staggered
            phivUU=v[i][j+1][k];        VOLUU=0.5*(ly[j+1]+ly[j]);  				//staggered
            phivUUU=v[i][j+2][k];       VOLUUU=0.5*(ly[j+2]+ly[j+1]);   			//staggered
            phivD=v[i][j-1][k];         VOLD=0.5*(ly[j-1]+ly[j-2]); 				//staggered
            phivDD=v[i][j-2][k];        VOLDD=0.5*(ly[j-2]+ly[j-3]);    			//staggered
            vSCHEMES();
        }
        }
    vsa=vadvel;
    }
    /*...Z-advection of V...*/
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            vadvel=vt;
        }
    else
        {
        if (wt>0)
        {
            phivU=v[i][j][k];       	VOLU=lz[k]; 						//collocated
            phivUU=v[i][j][k-1];        VOLUU=lz[k-1];  					//collocated
            phivUUU=v[i][j][k-2];       VOLUUU=lz[k-2]; 					//collocated
            phivD=v[i][j][k+1];         VOLD=lz[k+1];   					//collocated
            phivDD=v[i][j][k+2];        VOLDD=lz[k+2];  					//collocated
            vSCHEMES();
        }
        else
        {
            phivU=v[i][j][k+1];         VOLU=lz[k+1];   					//collocated
            phivUU=v[i][j][k+2];        VOLUU=lz[k+2];  					//collocated
            phivUUU=v[i][j][k+3];       VOLUUU=lz[k+3]; 					//collocated
            phivD=v[i][j][k];       	VOLD=lz[k]; 						//collocated
            phivDD=v[i][j][k-1];        VOLDD=lz[k-1];  					//collocated
            vSCHEMES();
        }
        }
    vta=vadvel;
    }
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            vadvel=vb;
        }
    else
        {
        if (wb>0)
        {
            phivU=v[i][j][k-1];         VOLU=lz[k-1];   					//collocated
            phivUU=v[i][j][k-2];        VOLUU=lz[k-2];  					//collocated
            phivUUU=v[i][j][k-3];       VOLUUU=lz[k-3]; 					//collocated
            phivD=v[i][j][k];       	VOLD=lz[k]; 						//collocated
            phivDD=v[i][j][k+1];        VOLDD=lz[k+1];  					//collocated
            vSCHEMES();
        }
        else
        {
            phivU=v[i][j][k];       	VOLU=lz[k]; 						//collocated
            phivUU=v[i][j][k+1];        VOLUU=lz[k+1];  					//collocated
            phivUUU=v[i][j][k+2];       VOLUUU=lz[k+2]; 					//collocated
            phivD=v[i][j][k-1];         VOLD=lz[k-1];   					//collocated
            phivDD=v[i][j][k-2];        VOLDD=lz[k-2];  					//collocated
            vSCHEMES();
        }
        }
    vba=vadvel;
    }
}

/*...FUNCTIONS NECESSARY FOR W-ADVECTION ROUTINE...*/
void wFOU()											//FIRST ORDER UPWIND (FOU)
{
    wadvel=phiwU;
}

void wSOU()											//SECOND ORDER UPWIND (SOU)
{
    wadvel=((2.0*VOLU+VOLUU)/(VOLU+VOLUU))*phiwU-(VOLU/(VOLU+VOLUU))*phiwUU;
}

void wQUICK()											//QUADRATIC UPSTREAM INTERPOLATION OF CONVECTION KINETICS (QUICK)
{
    wadvel=((VOLU*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLD+2.0*VOLU+VOLUU)))*phiwD+((VOLD*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLU+VOLUU)))*phiwU-((VOLD*VOLU)/((VOLU+VOLUU)*(VOLD+2.0*VOLU+VOLUU)))*phiwUU;
}

void wTVD()											//TOTAL VARIATION DIMINISHING (TVD)
{
    double xi=0.0;
    {
    xi=fabs((phiwU-phiwUU)/(phiwD-phiwUU));
        if (xi==0.0||xi==1.0)
        wadvel=phiwU;
        else if (xi>0.0&&xi<=0.3)
        wadvel=((3.0*VOLU+VOLUU)/(VOLU+VOLUU))*phiwU-((2.0*VOLU)/(VOLU+VOLUU))*phiwUU;
        else if (xi>0.3&&xi<=0.833)
        wadvel=((VOLU*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLD+2.0*VOLU+VOLUU)))*phiwD+((VOLD*(2.0*VOLU+VOLUU))/((VOLD+VOLU)*(VOLU+VOLUU)))*phiwU-((VOLD*VOLU)/((VOLU+VOLUU)*(VOLD+2.0*VOLU+VOLUU)))*phiwUU;
        else if (xi>0.833&&xi<1.0)
        wadvel=phiwD;
        else
        wadvel=phiwU;
    }
}

void wFiOU()											//FIFTH ORDER UPWIND (FiOU)
{
    //cell-center distances from the face
    double xU1=-0.5*VOLU;									//upwind cell-center
    double xU2=-(VOLU+0.5*VOLUU);								//up-upwind cell-center
    double xU3=-(VOLU+VOLUU+0.5*VOLUUU);							//up-up-upwind cell-center
    double xD1=0.5*VOLD;									//downwind cell-center
    double xD2=VOLD+0.5*VOLDD;									//down-downwind cell-center

    wadvel=(1.0-((xU1*xU3*xD1*xD2)/((xU2-xU1)*(xU2-xU3)*(xU2-xD1)*(xU2-xD2)))-((xU1*xU2*xD1*xD2)/((xU3-xU1)*(xU3-xU2)*(xU3-xD1)*(xU3-xD2)))-((xU1*xU2*xU3*xD2)/((xD1-xU1)*(xD1-xU2)*(xD1-xU3)*(xD1-xD2)))-((xU1*xU2*xU3*xD1)/((xD2-xU1)*(xD2-xU2)*(xD2-xU3)*(xD2-xD1))))*phiwU+((xU1*xU3*xD1*xD2)/((xU2-xU1)*(xU2-xU3)*(xU2-xD1)*(xU2-xD2)))*phiwUU+((xU1*xU2*xD1*xD2)/((xU3-xU1)*(xU3-xU2)*(xU3-xD1)*(xU3-xD2)))*phiwUUU+((xU1*xU2*xU3*xD2)/((xD1-xU1)*(xD1-xU2)*(xD1-xU3)*(xD1-xD2)))*phiwD+((xU1*xU2*xU3*xD1)/((xD2-xU1)*(xD2-xU2)*(xD2-xU3)*(xD2-xD1)))*phiwDD;
}

void wSCHEMES()											//sch=4 is CD2
{
    if (sch==1)
    wFOU();
    else if (sch==2)
    wSOU();
    else if (sch==3)
    wQUICK();
    else if (sch==5)
    wTVD();
    else
    wFiOU();
}

/*...ADVECTION ROUTINE FOR W-VELOCITY...*/
void wadrou()
{
    /*...X-advection of W...*/
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            wadvel=we;
        }
    else
        {
        if (ue>0)
        {
            phiwU=w[i][j][k];       	VOLU=lx[i]; 						//collocated
            phiwUU=w[i-1][j][k];        VOLUU=lx[i-1];  					//collocated
            phiwUUU=w[i-2][j][k];       VOLUUU=lx[i-2]; 					//collocated
            phiwD=w[i+1][j][k];         VOLD=lx[i+1];   					//collocated
            phiwDD=w[i+2][j][k];        VOLDD=lx[i+2];  					//collocated
            wSCHEMES();
        }
        else
        {
            phiwU=w[i+1][j][k];         VOLU=lx[i+1];   					//collocated
            phiwUU=w[i+2][j][k];        VOLUU=lx[i+2];  					//collocated
            phiwUUU=w[i+3][j][k];       VOLUUU=lx[i+3]; 					//collocated
            phiwD=w[i][j][k];       	VOLD=lx[i]; 						//collocated
            phiwDD=w[i-1][j][k];        VOLDD=lx[i-1];  					//collocated
            wSCHEMES();
        }
        }
    wea=wadvel;
    }
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            wadvel=ww;
        }
    else
        {
        if (uw>0)
        {
            phiwU=w[i-1][j][k];         VOLU=lx[i-1];   					//collocated
            phiwUU=w[i-2][j][k];        VOLUU=lx[i-2];  					//collocated
            phiwUUU=w[i-3][j][k];       VOLUUU=lx[i-3]; 					//collocated
            phiwD=w[i][j][k];       	VOLD=lx[i]; 						//collocated
            phiwDD=w[i+1][j][k];        VOLDD=lx[i+1];  					//collocated
            wSCHEMES();
        }
        else
        {
            phiwU=w[i][j][k];       	VOLU=lx[i]; 						//collocated
            phiwUU=w[i+1][j][k];        VOLUU=lx[i+1];  					//collocated
            phiwUUU=w[i+2][j][k];       VOLUUU=lx[i+2]; 					//collocated
            phiwD=w[i-1][j][k];         VOLD=lx[i-1];   					//collocated
            phiwDD=w[i-2][j][k];        VOLDD=lx[i-2];  					//collocated
            wSCHEMES();
        }
        }
    wwa=wadvel;
    }
    /*...Y-advection of W...*/
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            wadvel=wn;
        }
    else
        {
        if (vn>0)
        {
            phiwU=w[i][j][k];       	VOLU=ly[j]; 						//collocated
            phiwUU=w[i][j-1][k];        VOLUU=ly[j-1];  					//collocated
            phiwUUU=w[i][j-2][k];       VOLUUU=ly[j-2]; 					//collocated
            phiwD=w[i][j+1][k];         VOLD=ly[j+1];   					//collocated
            phiwDD=w[i][j+2][k];        VOLDD=ly[j+2];  					//collocated
            wSCHEMES();
        }
        else
        {
            phiwU=w[i][j+1][k];         VOLU=ly[j+1];   					//collocated
            phiwUU=w[i][j+2][k];        VOLUU=ly[j+2];  					//collocated
            phiwUUU=w[i][j+3][k];       VOLUUU=ly[j+3]; 					//collocated
            phiwD=w[i][j][k];       	VOLD=ly[j]; 						//collocated
            phiwDD=w[i][j-1][k];        VOLDD=ly[j-1];  					//collocated
            wSCHEMES();
        }
        }
    wna=wadvel;
    }
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            wadvel=ws;
        }
    else
        {
        if (vs>0)
        {
            phiwU=w[i][j-1][k];         VOLU=ly[j-1];   					//collocated
            phiwUU=w[i][j-2][k];        VOLUU=ly[j-2];  					//collocated
            phiwUUU=w[i][j-3][k];       VOLUUU=ly[j-3]; 					//collocated
            phiwD=w[i][j][k];       	VOLD=ly[j]; 						//collocated
            phiwDD=w[i][j+1][k];        VOLDD=ly[j+1];  					//collocated
            wSCHEMES();
        }
        else
        {
            phiwU=w[i][j][k];      	VOLU=ly[j]; 						//collocated
            phiwUU=w[i][j+1][k];        VOLUU=ly[j+1];  					//collocated
            phiwUUU=w[i][j+2][k];       VOLUUU=ly[j+2]; 					//collocated
            phiwD=w[i][j-1][k];         VOLD=ly[j-1];   					//collocated
            phiwDD=w[i][j-2][k];        VOLDD=ly[j-2];  					//collocated
            wSCHEMES();
        }
        }
    wsa=wadvel;
    }
    /*...Z-advection of W...*/
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            wadvel=wt;
        }
    else
        {
        if (wt>0)
        {
            phiwU=w[i][j][k];       	VOLU=0.5*(lz[k]+lz[k-1]);   				//staggered
            phiwUU=w[i][j][k-1];        VOLUU=0.5*(lz[k-1]+lz[k-2]);    			//staggered
            phiwUUU=w[i][j][k-2];       VOLUUU=0.5*(lz[k-2]+lz[k-3]);   			//staggered
            phiwD=w[i][j][k+1];         VOLD=0.5*(lz[k+1]+lz[k]);   				//staggered
            phiwDD=w[i][j][k+2];        VOLDD=0.5*(lz[k+2]+lz[k+1]);    			//staggered
            wSCHEMES();
        }
        else
        {
            phiwU=w[i][j][k+1];         VOLU=0.5*(lz[k+1]+lz[k]);   				//staggered
            phiwUU=w[i][j][k+2];        VOLUU=0.5*(lz[k+2]+lz[k+1]);    			//staggered
            phiwUUU=w[i][j][k+3];       VOLUUU=0.5*(lz[k+3]+lz[k+2]);   			//staggered
            phiwD=w[i][j][k];       	VOLD=0.5*(lz[k]+lz[k-1]);   				//staggered
            phiwDD=w[i][j][k-1];        VOLDD=0.5*(lz[k-1]+lz[k-2]);    			//staggered
            wSCHEMES();
        }
        }
    wta=wadvel;
    }
    {
    if (sch==4)											//CENTRAL DIFFERENCE
        {
            wadvel=wb;
        }
    else
        {
        if (wb>0)
        {
            phiwU=w[i][j][k-1];         VOLU=0.5*(lz[k-1]+lz[k-2]); 				//staggered
            phiwUU=w[i][j][k-2];        VOLUU=0.5*(lz[k-2]+lz[k-3]);    			//staggered
            phiwUUU=w[i][j][k-3];       VOLUUU=0.5*(lz[k-3]+lz[k-4]);   			//staggered
            phiwD=w[i][j][k];       	VOLD=0.5*(lz[k]+lz[k-1]);   				//staggered
            phiwDD=w[i][j][k+1];        VOLDD=0.5*(lz[k+1]+lz[k]);  				//staggered
            wSCHEMES();
        }
        else
        {
            phiwU=w[i][j][k];       	VOLU=0.5*(lz[k]+lz[k-1]);   				//staggered
            phiwUU=w[i][j][k+1];        VOLUU=0.5*(lz[k+1]+lz[k]);  				//staggered
            phiwUUU=w[i][j][k+2];       VOLUUU=0.5*(lz[k+2]+lz[k+1]);   			//staggered
            phiwD=w[i][j][k-1];         VOLD=0.5*(lz[k-1]+lz[k-2]); 				//staggered
            phiwDD=w[i][j][k-2];        VOLDD=0.5*(lz[k-2]+lz[k-3]);    			//staggered
            wSCHEMES();
        }
        }
    wba=wadvel;
    }
}

void mustar_FX()
{
        mue=must[i][j][k];									//stagger
        muw=must[i-1][j][k];									//stagger
            double mu_Ucv_P=(must[i][j][k]*must[i-1][j][k])/(((0.5*lx[i-1]*must[i][j][k])+(0.5*lx[i]*must[i-1][j][k]))/(0.5*(lx[i-1]+lx[i])));
        {
            double mu_Ucv_N=(must[i][j+1][k]*must[i-1][j+1][k])/(((0.5*lx[i-1]*must[i][j+1][k])+(0.5*lx[i]*must[i-1][j+1][k]))/(0.5*(lx[i-1]+lx[i])));
            double mu_Ucv_S=(must[i][j-1][k]*must[i-1][j-1][k])/(((0.5*lx[i-1]*must[i][j-1][k])+(0.5*lx[i]*must[i-1][j-1][k]))/(0.5*(lx[i-1]+lx[i])));
            mun=(mu_Ucv_N*mu_Ucv_P)/(((0.5*ly[j]*mu_Ucv_N)+(0.5*ly[j+1]*mu_Ucv_P))/(0.5*(ly[j]+ly[j+1])));
            mus=(mu_Ucv_P*mu_Ucv_S)/(((0.5*ly[j-1]*mu_Ucv_P)+(0.5*ly[j]*mu_Ucv_S))/(0.5*(ly[j-1]+ly[j])));
        }
        {
            double mu_Ucv_T=(must[i][j][k+1]*must[i-1][j][k+1])/(((0.5*lx[i-1]*must[i][j][k+1])+(0.5*lx[i]*must[i-1][j][k+1]))/(0.5*(lx[i-1]+lx[i])));
            double mu_Ucv_B=(must[i][j][k-1]*must[i-1][j][k-1])/(((0.5*lx[i-1]*must[i][j][k-1])+(0.5*lx[i]*must[i-1][j][k-1]))/(0.5*(lx[i-1]+lx[i])));
            mut=(mu_Ucv_T*mu_Ucv_P)/(((0.5*lz[k]*mu_Ucv_T)+(0.5*lz[k+1]*mu_Ucv_P))/(0.5*(lz[k]+lz[k+1])));
            mub=(mu_Ucv_P*mu_Ucv_B)/(((0.5*lz[k-1]*mu_Ucv_P)+(0.5*lz[k]*mu_Ucv_B))/(0.5*(lz[k-1]+lz[k])));
        }
}

void mustar_FY()
{
            double mu_Vcv_P=(must[i][j][k]*must[i][j-1][k])/(((0.5*ly[j-1]*must[i][j][k])+(0.5*ly[j]*must[i][j-1][k]))/(0.5*(ly[j-1]+ly[j])));
        {
            double mu_Vcv_E=(must[i+1][j][k]*must[i+1][j-1][k])/(((0.5*ly[j-1]*must[i+1][j][k])+(0.5*ly[j]*must[i+1][j-1][k]))/(0.5*(ly[j-1]+ly[j])));
            double mu_Vcv_W=(must[i-1][j][k]*must[i-1][j-1][k])/(((0.5*ly[j-1]*must[i-1][j][k])+(0.5*ly[j]*must[i-1][j-1][k]))/(0.5*(ly[j-1]+ly[j])));
            mue=(mu_Vcv_E*mu_Vcv_P)/(((0.5*lx[i]*mu_Vcv_E)+(0.5*lx[i+1]*mu_Vcv_P))/(0.5*(lx[i]+lx[i+1])));
            muw=(mu_Vcv_P*mu_Vcv_W)/(((0.5*lx[i-1]*mu_Vcv_P)+(0.5*lx[i]*mu_Vcv_W))/(0.5*(lx[i-1]+lx[i])));
        }
        mun=must[i][j][k];									//stagger
        mus=must[i][j-1][k];									//stagger
        {
            double mu_Vcv_T=(must[i][j][k+1]*must[i][j-1][k+1])/(((0.5*ly[j-1]*must[i][j][k+1])+(0.5*ly[j]*must[i][j-1][k+1]))/(0.5*(ly[j-1]+ly[j])));
            double mu_Vcv_B=(must[i][j][k-1]*must[i][j-1][k-1])/(((0.5*ly[j-1]*must[i][j][k-1])+(0.5*ly[j]*must[i][j-1][k-1]))/(0.5*(ly[j-1]+ly[j])));
            mut=(mu_Vcv_T*mu_Vcv_P)/(((0.5*lz[k]*mu_Vcv_T)+(0.5*lz[k+1]*mu_Vcv_P))/(0.5*(lz[k]+lz[k+1])));
            mub=(mu_Vcv_P*mu_Vcv_B)/(((0.5*lz[k-1]*mu_Vcv_P)+(0.5*lz[k]*mu_Vcv_B))/(0.5*(lz[k-1]+lz[k])));
        }
}

void mustar_FZ()
{
            double mu_Wcv_P=(must[i][j][k]*must[i][j][k-1])/(((0.5*lz[k-1]*must[i][j][k])+(0.5*lz[k]*must[i][j][k-1]))/(0.5*(lz[k-1]+lz[k])));
        {
            double mu_Wcv_E=(must[i+1][j][k]*must[i+1][j][k-1])/(((0.5*lz[k-1]*must[i+1][j][k])+(0.5*lz[k]*must[i+1][j][k-1]))/(0.5*(lz[k-1]+lz[k])));
            double mu_Wcv_W=(must[i-1][j][k]*must[i-1][j][k-1])/(((0.5*lz[k-1]*must[i-1][j][k])+(0.5*lz[k]*must[i-1][j][k-1]))/(0.5*(lz[k-1]+lz[k])));
            mue=(mu_Wcv_E*mu_Wcv_P)/(((0.5*lx[i]*mu_Wcv_E)+(0.5*lx[i+1]*mu_Wcv_P))/(0.5*(lx[i]+lx[i+1])));
            muw=(mu_Wcv_P*mu_Wcv_W)/(((0.5*lx[i-1]*mu_Wcv_P)+(0.5*lx[i]*mu_Wcv_W))/(0.5*(lx[i-1]+lx[i])));
        }
        {
            double mu_Wcv_N=(must[i][j+1][k]*must[i][j+1][k-1])/(((0.5*lz[k-1]*must[i][j+1][k])+(0.5*lz[k]*must[i][j+1][k-1]))/(0.5*(lz[k-1]+lz[k])));
            double mu_Wcv_S=(must[i][j-1][k]*must[i][j-1][k-1])/(((0.5*lz[k-1]*must[i][j-1][k])+(0.5*lz[k]*must[i][j-1][k-1]))/(0.5*(lz[k-1]+lz[k])));
            mun=(mu_Wcv_N*mu_Wcv_P)/(((0.5*ly[j]*mu_Wcv_N)+(0.5*ly[j+1]*mu_Wcv_P))/(0.5*(ly[j]+ly[j+1])));
            mus=(mu_Wcv_P*mu_Wcv_S)/(((0.5*ly[j-1]*mu_Wcv_P)+(0.5*ly[j]*mu_Wcv_S))/(0.5*(ly[j-1]+ly[j])));
        }
        mut=must[i][j][k];									//stagger
        mub=must[i][j][k-1];									//stagger
}
