#include "main.h"

#define LightSpeed 3e8
#define PI 3.1415926
#define MAX(x,y) ((x>y)?x:y)

#define Pi 3.1415926
#define RAD2DEG (180./Pi)
#define DEG2RAD (Pi/180.)


double get_expected_phase(double Az,double El,double Lambda,double D,double H)
 {
  double res;
/// calc expected phase for given Azimuth and Elevation
  res=2.*3.1415*1.0/Lambda*(D*sqrt(cos(Az)*cos(Az)-sin(El)*sin(El))+H*sin(El));
  return res;

//  return atan2(sin(res/2.0)/cos(res/2.0),1.0)*2.0;

 }

double phi_to_elev(double phi,double Az,double Lambda,double D,double H)
{
#define WRONG_RESULT -999.0
 double M=phi/(2.*3.1415)*Lambda/D;

 double det=(cos(Az)*cos(Az)-M*M)*(1+H*H/(D*D))+M*M*H*H/(D*D);
 if(det>0.)
  {
   double res=(sqrt(det)+H/D)/(1+H*H/(D*D));
   if (fabs(res)<=1)
    return asin(res);
   else
    return WRONG_RESULT;
  }
 return WRONG_RESULT;
}

// #include "read_config.h"


int main(int argc,char *argv[]) {
  FILE *fp;
  struct RadarParm prm;
  struct RadarParm *PRM;
  struct FitData fit;
  struct FitData *FIT;
  char* radarname;
  double Xlat0=(double)56.5,Xlong0=(double)58.5;
  float DayNo=0;
  long SELECT_BEAM=0;


//  config_record_type* config_record;
//  long config_record_len;

#ifndef WRONG_PHASE
#define WRONG_PHASE -999.9
#endif

  double phase=WRONG_PHASE;
//   read_config_phase(&config_record,&config_record_len);


  
  if(argc<2)
   {
      fprintf(stderr,"usage %s fname outfname\n",argv[0]);
      exit(1);
   }
  fp=fopen(argv[1],"r");
  long CHANNEL=-1;//atol(argv[2]);
  long BEAM=-1;//atol(argv[3]);
  FILE* fpout;
  if(argv[2])
   fpout=fopen(argv[2],"wb");
  
  if (fp==NULL) 
   {
    fprintf(stderr,"File %s not found.\n",argv[1]);
    exit(-1);
   }

  PRM=RadarParmMake();
  FIT=FitMake();
  prm=*PRM;
  fit=*FIT;

  int init=0;

double range;
double phi0;
double TdiffNew=0.0;
double PrmIntX=0.0; 

int front=1;
double alias=0.0;
double elv_out[1000];
double j=0;
  while(FitFread(fp,&prm,&fit) !=-1) 
   {
     long long int curdate;

     curdate=(long long int)prm.time.yr;
     curdate=curdate*100+(long long int)prm.time.mo;
     curdate=curdate*100+(long long int)prm.time.dy;
     curdate=curdate*100+(long long int)prm.time.hr;
     curdate=curdate*100+(long long int)prm.time.mt;
     
     if(fit.elv!=NULL)
       {
    
long i;

     for(i=0;i<prm.nrang;i++)
      {
       fit.elv[i].normal=
       fit.elv[i].high=
       fit.elv[i].low=0.;
       
       if(
        //berng fit.rng[i].p_l>3 && 
        fit.rng[i].qflg>0 && 
        fit.xrng[i].phi0_err>0 
//berng  && prm.channel==CHANNEL 
//       && (long)prm.bmnum==BEAM
       )
       {
    	double A,B,C,expected_phi;
    	double D;
    	double H;  ///RawACF reverse sign for interferometer, we must revers height from -3.8 to 3.8
//KER
	D=-89.6; H=0.0;
	
	double D1,beam,expected_elev;
	double Lambda;
	Lambda=150000./(double)prm.tfreq*2.0;

	beam=(double)prm.bmnum;

	A=7.9064+0.0012254*(double)prm.tfreq; 

    	B=0;
    	C=0;

	double BeamWidth=3.24;
	double Az;

	Az=(beam-7.5)*BeamWidth*3.1415/180.;
	D1=D*cos(Az);

	range=prm.frang+prm.rsep*(double)i;
	double Re=6371;
	double altitude=90.;
	expected_elev=asin((((Re+altitude)*(Re+altitude)-(Re*Re+range*range))/(2.*Re*range)));
	double expected_elev_plain;
	expected_elev_plain=atan2(altitude,sqrt(range*range-altitude*altitude));
	expected_phi=get_expected_phase(Az,expected_elev,Lambda,D,H);

	double NormD1Lambda;
	NormD1Lambda=D1/Lambda;//-floor(D1/Lambda); /// Fix according to Shepherd 2017
	
	
	double corrected_phi;
	double fr_norm;
	fr_norm=(double)prm.tfreq/10000.;

	corrected_phi=fit.xrng[i].phi0-A;//-B*fr_norm-C*fr_norm*fr_norm;  //for KER

	double calibrated_elev;
	double tmp_corr_phi;
	long n;
	double min_calibrated_elev=100;
	double min_tmp_corr_phi=100;
	long min_n=10;

	double max_calibrated_elev=-100;
	double max_tmp_corr_phi=-100;
	long max_n=-10;
	for(n=5;n>-5;n--)
	{
	 tmp_corr_phi=corrected_phi+(double)n*2.*3.1415;
	 calibrated_elev=phi_to_elev(tmp_corr_phi,Az,Lambda,D,H);

        if(
    	calibrated_elev > WRONG_RESULT
//    	&&
//	calibrated_elev*180./3.14 > 7.	

        && range>0. 
    	&& beam>=0 && beam<=15
        )
        {
         if(calibrated_elev<min_calibrated_elev)
          {
           max_calibrated_elev=min_calibrated_elev;
           max_tmp_corr_phi=min_tmp_corr_phi;
           max_n=min_n;
           min_calibrated_elev=calibrated_elev;
           min_tmp_corr_phi=tmp_corr_phi;
           min_n=n;
          }
        }
        }

	if(min_n<10)
	{

	    fit.elv[i].normal=min_calibrated_elev*180./3.14;
	    fit.elv[i].high=max_calibrated_elev*180./3.14;
    	    fit.elv[i].low=min_calibrated_elev*180./3.14;
    	}
	char substring[255];
	substring[0]=0;

#ifdef DEBUG
	if(min_n<10 
	 && (double)prm.frang+(double)prm.rsep*i<350
//	 && (long)prm.bmnum >=5 && (long)prm.bmnum <=10
//	&& prm.tfreq<12000
    		
	 )
         printf("%lf"
        	"\t%.0lf\t%.3lf"
        	"\t%ld %ld %ld"
        	"\t%.0lf %.0lf %.0lf"
        	"\t%ld"
        	"\t%lf %lf"
        	"\t%lf %lf"
        	"\t%.0lf %ld"
        	"\t%lf %lf %lf"
        	"\t%lf"
        	"\t%s"
        	"\t%f"
        	"\t%lf %lf %ld"
        	"\t\n",
			(double)prm.time.hr+(double)prm.time.mt/60.+(double)prm.time.sc/3600.,

    			(double)prm.frang+(double)prm.rsep*i,
    			min_calibrated_elev*180./3.14,

    			(long)prm.channel,
    			(long)prm.tfreq,
    			(long)prm.bmnum,

    			(double)fit.rng[i].p_l,
    			(double)fit.rng[i].v,
    			(double)fit.rng[i].w_l,

    			(long)fit.rng[i].gsct,
    			
    			expected_phi,
    			corrected_phi,

    			fit.elv[i].high,
    			fit.elv[i].low,
    			
    			sin(min_calibrated_elev)*(range+(double)prm.rsep/2.)+(range*range)/(2.*Re),
    			min_n,
    			
    			fit.xrng[i].phi0,
    			fit.rng[i].phi0,
    			phase,
    			
    			A,
    			substring,
			expected_elev*180./3.14,

			max_calibrated_elev*180./3.14,
    			sin(max_calibrated_elev)*(range+(double)prm.rsep/2.)+(range*range)/(2.*Re),
    			max_n

    			);
#endif

	}
      }
     }
    if(fpout)
     FitFwrite(fpout,&prm,&fit);
   }
 
  fclose(fp);
  if(fpout)
   fclose(fpout);
  

}

