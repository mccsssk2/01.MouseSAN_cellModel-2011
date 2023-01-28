// Ca subsystem analysis
// since my volumes are different to Maltsev, I need to 
// do this analysis, at least roughly. thank god I have NGS.

#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<malloc.h>
#include<time.h>

#define tracesornot 0 // 0 for no traces, 1 for traces.
#define pupornot 0 // 0 is pumpvmax, 1 is pup

/*
For Random Numbers
*/

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/****************************************************************************************/


#define ddt 0.001 // ms time step.

#define capacitance 0.025 // nF

#define vcell 3.0
#define vrel 0.0036
#define vsub 0.03328117

#define vup  0.0348
#define vi 1.34671883

#define Mgi 2.5	      		// Internal Mg2+ concentration (fixed)

/*
Calcium handling
*/

#define tdifca 0.04 // as in Maltsev. This is diffusion from Cai to Casub, this is not to do with the transients.

// #define Pup 0.025 // for a frequency of 5 hertz or more.

#define Kup 0.0006 	// 0.0006
#define nup 1.0	//% SR Ca2+ uptake Km (mM) and Hill coefficient
#define Ttr 40.0 	// Time constant for Ca2+ transfer from NSR to JSR (ms) 60

//------------------------------------------------------------------------------------------------------------------------------
#define ConcTC 0.031            //              % Concentration of Troponin-Ca complex (mM)
#define ConcTMC 0.062           //      % Concentration of Troponin-Mg complex (mM)
#define kfTC 88.8
#define kfTMC 237.7             //% Rate constant for Ca2+ binding to Troponin(mM/ms)
#define kbTC 0.446
#define kbTMC 0.00751   //      % Rate constant for Ca2+ unbinding from Troponin (ms-1)
#define kfTMM 2.277             //              % Rate constant for Mg2+ binding to Troponin(mM/ms)
#define kbTMM 0.751             //              % Rate constant for Mg2+ unbinding from Troponin (ms-1)
//------------------------------------------------------------------------------------------------------------------------------
#define ConcCM 0.045    //                      % Concentration of Calmodulin (mM)
#define kfCM 227.7              // manip.               % Rate constant for Ca2+ binding to Calmodulin (mM/ms)
#define kbCM 0.542              //              % Rate constant for Ca2+ unbinding from Calmodulin (mM/ms)
//------------------------------------------------------------------------------------------------------------------------------
#define ConcCQ 10.0             // % Concentration of Calsequestrin (mM)
#define kfCQ 0.534              //              % Rate constant for Ca2+ binding to Calsequestrin (mM/ms)
#define kbCQ 0.445              //              % Rate constant for Ca2+ unbinding from Calsequestrin (mM/ms)
//===========================================================================


/*
RyR Markov chain
*/

#define koca 10.0  // mM-2 ms-1
#define kom 0.06 // ms-1
#define kica 0.5 // mM-1ms-1
#define kim 0.005 // ms-1
#define eca50sr 0.45 // mM
// #define ks 250000.0 // ms-1
#define maxsr 15.0
#define minsr 1.0
#define hsrr 2.5

/*
The new pump
 */

#define pumpkmf 0.000246 // mM
#define pumpkmr 1.7 // mM
#define pumphill 1.787


/*
Ionic homeostasis
*/

double ca_flux;

/*
Calcium diffusion
*/

double Jcadif;

double Jrel, Jup, Jtr,carel,caup;
double cai;             // mM
double casub;           // mM 


double dFtc,Ftc;
double dFtmc,Ftmc,Ftmm;
double dFtmm,Ftmm;
double dFcms,Fcms;
double dFcmi,Fcmi;
double dFcq,Fcq;

/*
RyR markov chain
*/

double resting, open, inactivated, resting_inactivated;
double kcasr,kosrca,kisrca;

/*
the SERCA pump
*/

double pumpvmax = 0.001; // mM/ms: replacement of Pup

double sstime;

double total_time = 0.0;

double carelmin[1000], carelmax[1000], caupmin[1000], caupmax[1000], caimin[1000], caimax[1000], casubmin[1000], casubmax[1000];
double casubmintime[1000];

int counter = 0, measureCounter = 0, measureCountermax = 0;

double dvdt = 0.0, dvdtold = 0.0;

int filecounter = 0;
int oscillatingOrnot = 0;

double Pup, ks;

int solution_broke = 0, somerandomnumber = 23334556;

/*
The random number generator.
*/

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

int main(){

FILE *outputcurrents,*parameters, *random;
char *str;

/*
Initialising the random number generator is important.
*/

// init_by_array(init, length);

random = fopen("/dev/urandom","rb");
fread(&somerandomnumber, sizeof(unsigned int), 1, random);
fclose(random);

// somerandomnumber = rand();
init_genrand(somerandomnumber);

for(filecounter = 0; filecounter < 50000; filecounter++){

Pup = genrand_real1()*0.1;
pumpvmax = genrand_real1()*0.1; // around 0.0001 in the Shannon model.
ks = genrand_real1()*1800000.0;

// Pup = 0.04;
// ks = 1300000.0;



 if(tracesornot==1){
 str = malloc(32*sizeof(char));
 sprintf(str,"central%d.dat",filecounter);
 outputcurrents = fopen(str,"w");
 free(str);
 }

// initial conditions

carel = 0.029;
caup = 1.35;
casub = 0.000223;
Ftc = 0.02;
Ftmc = 0.22;
Ftmm = 0.69;
Fcms = 0.089;
Fcmi = 0.042;
Fcq = 0.032;
cai = 0.0001;

resting = 0.7499;
open = 0.0000034;
inactivated = 0.0000011;
resting_inactivated = 0.25;

ca_flux = 0.0; // no ionic currents

total_time = 12000.0;

counter = 0;

for(measureCounter = 0; measureCounter < 1000; measureCounter++){
carelmin[measureCounter] = 10000.0;
carelmax[measureCounter] = -10000.0;
caupmin[measureCounter] = 10000.0;
caupmax[measureCounter] = -10000.0;
caimin[measureCounter] = 10000.0;
caimax[measureCounter] = -10000.0;
casubmin[measureCounter] = 10000.0;
casubmax[measureCounter] = -10000.0;
casubmintime[measureCounter] = -10000.0;
}

measureCounter = 0;

dvdt = 0.0; dvdtold = 0.0;

for(sstime=0.0;sstime<=total_time;sstime = sstime + ddt){         // time loop starts here

counter++;

if(cai!=cai){
printf("solution broke %f\n",sstime);
 break;
 }

/*****************************************************************************/

/* Ca2+ difusion */
Jcadif = (casub - cai)/tdifca;

/* Ca handling in SR (mM/ms) */
   
 kcasr = maxsr - (maxsr - minsr)/(1.0 + pow(eca50sr/carel,hsrr));
 kosrca = koca/kcasr;
 kisrca = kica*kcasr;
 resting = resting + ddt*(kim*resting_inactivated - kisrca*casub*resting 
 - kosrca*casub*casub*resting + kom*open);
 open = open + ddt*(kosrca*casub*casub*resting 
 - kom*open - kisrca*casub*open + kim*inactivated);
 inactivated = inactivated + ddt*(kisrca*casub*open 
 - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated);
 resting_inactivated = resting_inactivated + ddt*(kom*inactivated 
 - kosrca*casub*casub*resting_inactivated - kim*resting_inactivated + kisrca*casub*resting);

Jrel = ks*open*(carel - casub);
    
 if(pupornot==1)
    Jup  = Pup/(1.0+Kup/cai);
 else
Jup = pumpvmax*(pow(cai/pumpkmf,pumphill) - pow(casub/pumpkmr,pumphill))/(1.0 + pow(cai/pumpkmf,pumphill) + pow(casub/pumpkmr,pumphill));
            
    Jtr  = (caup - carel)/Ttr;

    dFtc  = kfTC*cai*(1.0-Ftc)-kbTC*Ftc;
    dFtmc = kfTMC*cai*(1.0-Ftmc-Ftmm)-kbTMC*Ftmc;
    dFtmm = kfTMM*Mgi*(1.0-Ftmc-Ftmm)-kbTMM*Ftmm;

    dFcms = kfCM*casub*(1.0-Fcms)-kbCM*Fcms;
    dFcmi = kfCM*cai*(1.0-Fcmi)-kbCM*Fcmi;
    dFcq  = kfCQ*carel*(1.0-Fcq)-kbCQ*Fcq;

dvdtold = dvdt;

dvdt = ((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);

 casub = casub + ddt*((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);
 cai = cai + ddt*((Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmi + ConcTC*dFtc + ConcTMC*dFtmc));

 carel = carel + ddt*(Jtr - Jrel - ConcCQ*dFcq);
 caup = caup + ddt*(Jup-Jtr*vrel/vup);

    Ftc = Ftc + ddt*dFtc;                   //	    % dfTnCa/dt
    Ftmc = Ftmc + ddt*dFtmc;                  //	% dfTnMgCa/dt
    Ftmm = Ftmm + ddt*dFtmm;                  //	% dfTnMgMg/dt
    Fcms = Fcms + ddt*dFcms;                 //	    % dfCMs/dt
    Fcmi = Fcmi + ddt*dFcmi;                 //	    % dfCMi/dt
    Fcq = Fcq + ddt*dFcq;                    //	    % dfCQ/dt

/*************************************************************************************************/

// measurement

if(sstime>=total_time*0.6){

if(casubmin[measureCounter]>0&&dvdtold<0&&dvdt>=0){
measureCounter++;
 casubmin[measureCounter] = casub; // this is the first minimum, so count from here.
 casubmintime[measureCounter] = sstime;
 }

if(casubmax[measureCounter]<casub){
 casubmax[measureCounter] = casub;
 }

if(carelmin[measureCounter]>carel)
carelmin[measureCounter]=carel;

if(carelmax[measureCounter]<carel)
carelmax[measureCounter]=carel;

if(caupmin[measureCounter]>caup)
caupmin[measureCounter]=caup;

if(caupmax[measureCounter]<caup)
caupmax[measureCounter]=caup;

if(caimin[measureCounter]>cai)
caimin[measureCounter]=cai;

if(caimax[measureCounter]<cai)
caimax[measureCounter]=cai;

} // end of measurement

if(measureCounter>=998) break;

 if(tracesornot==1)
   if(counter>1000&&sstime>(total_time-2000)){

fprintf(outputcurrents,"%20.10f\t",sstime); 
fprintf(outputcurrents,"%20.10f\t",carel); 
fprintf(outputcurrents,"%20.10f\t",caup); 
fprintf(outputcurrents,"%20.10f\t",casub); 
fprintf(outputcurrents,"%20.10f\t",cai); 
fprintf(outputcurrents,"%20.10f\t",Ftc); 
fprintf(outputcurrents,"%20.10f\t",Ftmc); 
fprintf(outputcurrents,"%20.10f\t",Ftmm); 
fprintf(outputcurrents,"%20.10f\t",Fcms); 
fprintf(outputcurrents,"%20.10f\t",Fcmi); 
fprintf(outputcurrents,"%20.10f\t",Fcq); 
fprintf(outputcurrents,"%20.10f\t",resting); 
fprintf(outputcurrents,"%20.10f\t",open);
fprintf(outputcurrents,"%20.10f\t",inactivated); 
fprintf(outputcurrents,"%20.10f\t",resting_inactivated); 
fprintf(outputcurrents,"%20.10f\t",Jrel);
fprintf(outputcurrents,"%20.10f\t",Jup); 
fprintf(outputcurrents,"%20.10f\t",Jtr); 

fprintf(outputcurrents,"\n");

/*
printf("%20.10f\t",sstime); 
printf("%20.10f\t",carel); 
printf("%20.10f\t",caup); 
printf("%20.10f\t",casub); 
printf("%20.10f\t",cai); 
printf("%20.10f\t",Ftc); 
printf("%20.10f\t",Ftmc); 
printf("%20.10f\t",Ftmm); 
printf("%20.10f\t",Fcms); 
printf("%20.10f\t",Fcmi); 
printf("%20.10f\t",Fcq); 
printf("%20.10f\t",resting); 
printf("%20.10f\t",open);
printf("%20.10f\t",inactivated); 
printf("%20.10f\t",resting_inactivated); 


printf("\n");

*/


counter = 0;
}

} // end of time loop.

 if(tracesornot==1)
fclose(outputcurrents);

measureCountermax = measureCounter;
oscillatingOrnot = 1; // 1 if oscillating

if(measureCountermax==0){

oscillatingOrnot = 0;

measureCountermax = 10;
measureCounter = measureCountermax-1;

casubmintime[measureCounter] = -20; 

casubmin[measureCounter] = casub; 
casubmax[measureCounter] = casub;

carelmin[measureCounter]=carel;
carelmax[measureCounter]=carel;

caupmin[measureCounter]=caup;
caupmax[measureCounter]=caup;

caimin[measureCounter]=cai;
caimax[measureCounter]=cai;
 
}

 str = malloc(32*sizeof(char));
 sprintf(str,"Params.dat");
 parameters = fopen(str,"a");
 free(str);



for(measureCounter = measureCountermax-1; measureCounter < measureCountermax; measureCounter++){

  if(fabs(casubmintime[measureCounter]-casubmintime[measureCounter-1])>1500.0){
    casubmintime[measureCounter] = 0.0;
    casubmintime[measureCounter-1] = 0.0;
  }

fprintf(parameters,"%d\t%d\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\t%20.10f\n",filecounter,measureCounter,pumpvmax,ks,Pup,
casubmintime[measureCounter]-casubmintime[measureCounter-1],
casubmax[measureCounter]-casubmin[measureCounter],
	caimax[measureCounter]-caimin[measureCounter]);
 }

fclose(parameters);

} // end of fileloopcounter.

return 0;
}
