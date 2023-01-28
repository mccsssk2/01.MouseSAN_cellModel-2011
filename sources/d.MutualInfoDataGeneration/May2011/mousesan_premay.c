/*
This is the one on 25 Feb. 2011.
To do:
IKr TC.
INa1.1 more important.
ICaL1.2 more important.
Ks and Pup in good range - I have too much.
ISO phase in as simulation.
*/

#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<malloc.h>
#include<time.h>

#define if_ap_profiles 1 // 0 = no profiles, 1 = 2 sec. profiles.
#define screen_output 1 // 0 = no screen output, 1 yes screen output.

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

#define R 8.314472 // units:
#define T 310.5    // Kelvins
#define F 96.4846  // units:

#define number_of_apds 5000

/*
THIS IS OLD:
0 = 100% block only simulation.
1 = gh from 0% to 500%.
2 = gNattxr from 0% to 500%.
3 = gcal12
4 = gcal13
5 = gcat
6 = gst
NO 7 to 9
10 = gna_ttxs
11 = gkr
NO 12 or 13
14 = kNaCa
15 = Prel
17 = Pup
18 = If-v/12
*/

#define whichParamForBif 0

/***************************************************************************************/

#define ddt 0.001 // ms time step.

#define capacitance 0.025 // nF

// the volumes are now correct with proportions and differences.
#define vcell 3.0 // pL

#define l_cell 66.3767257 // microM
#define r_cell 3.792956 // microM

#define vrel 0.0036
#define vsub 0.03328117
#define vup  0.0348
#define vi 1.34671883

#define Mgi 2.5	 	// Internal Mg2+ concentration
#define nao 140.0       // mM
#define cao 1.8         // mM
#define ko 5.4          // mM

/*
Ist parameters 
*/

#define gst 0.00006 // it can be lower. Free parameter. Maltsev has 0.000075 microS
#define eist 17.0 // mV reversal of Ist

/*
Ib parameters
*/

#define gbna 0.0001215         //uS 
#define gbca 0.000015         //uS 
#define gbk  0.0000025         //uS

/*
IK1 parameters
*/

#define	gk1 0.229*0.0039228*0.9             // microS: mangoni is 0.0009 microS

/*
IKr parameters
*/

// #define gkr 2.2*0.002955      // micro-S   increased it temporarily to improve AP

/*
IKs parameters
*/

#define gks 0.000299    // micro-S ki=140 mM

/*
ICaL 1.2 and 1.3 parameters.
*/

// the standard multiplication factor is 4. Anything above that, and I need to increase the ICaT also.
/* 
#define icals_scaling 4.0
#define gcal12 0.0010*4.0*1.75
#define gcal13 0.0030*4.0*2.5
*/
 
#define ecal 47.0
#define kmfca 0.00035
#define alpha_fca 0.021

/*
ICaT
*/

#define all_ica_multiplier 1.0
#define ecat 45.0

/*
INa Nav 1.1 (TTXS) and 1.5 (TTXR)
*/

#define enattxr  41.5761

/*
If parameters
*/

#define multiplier2 1.0

/*
The conductance of gh from IV is 0.009 micro S. However, Mangoni has gh = 0.00324.
This is more relevant during AP, as we have stronger calcium handling.
gh can be increased by a factor of 2 if required.
*/

/*
Isus parameters
*/

#define gsus 0.00039060        // micro-S

/*
Ito parameters
*/

/*
INaK parameters
*/

// the higher f is, the lower Nai becomes eventually. to increase Nai, reduce pump, to reduce Nai, increase pump.

// Kurata/Maltsev INaK model
#define inakmax_multiplier 1.85
#define inakmax inakmax_multiplier*0.077 // 2.88 pA/pF as in Kurata and Maltsev x conversion to give nA
#define kmnap 14.0 // mM
#define kmkp 1.4 // mM

// % INaCa
// Kurata central value is 3.125 nA for 25 pF cell.
// #define kNaCa 5.5
#define K1ni 395.3
#define K1no 1628
#define K2ni 2.289
#define K2no 561.4
#define K3ni 26.44
#define K3no 4.663
#define Kci 0.0207
#define Kco 3.663
#define Kcni 26.44
#define Qci 0.1369
#define Qco 0.0
#define Qn 0.4315

/*
Calcium handling
*/

#define tdifca 0.04 // as in Maltsev. This is diffusion from Cai to Casub, this is not to do with the transients.
#define Prel 2.5 	//	% SR Ca2+ release rate constant (ms-1): as given in the word document

/*
These values of Pup and Ks give a Ca2+ CL of 195 ms, and good profiles 
for Cai, Casub, Jrel etc. Now I need to keep this sane, and REDUCE Pup eventually to get into the
good "CL" range.
At Pup = 0.03 and Ks = 600000, its an oscillator.
*/

#define Krel 0.0015
#define nrel 2.0	//% SR Ca2+ release Km (mM) and Hill coefficient
#define Kup 0.0006 	// 0.0006
#define nup 1.0		//% SR Ca2+ uptake Km (mM) and Hill coefficient
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
#define kfCM 237.7              // manip.               % Rate constant for Ca2+ binding to Calmodulin (mM/ms)
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
// #define ks 600000.0 // ms-1 for a freuqcy of 5 hertz or more
#define maxsr 15.0
#define minsr 1.0
#define hsrr 2.5

/*
Model variables
*/

double v, vnew,sstime;
int  current_trace_counter = 0, output_counter = 0;

double vclamp;
double total_current;

/*
Ion concentrations and reversal potentials
*/

double cai;             // mM
double casub;           // mM 

double ena,eca,ek,eks;

/*
Ist variables.
*/

double ist;
double qa, qi, tauqa,tauqi,alphaqa,betaqa,alphaqi,betaqi;
double dst; // ist gating variables
double fst;

/*
Ib variables
*/

double ibca,ibna,ibk,ib;

/*
IK1 variables
*/

double ik1 ;
double xk1inf;

/*
ICaT variables
*/

double icat;
double dt;
double ft;
double dt_inf, tau_dt, ft_inf, tau_ft;

/*
IKr variables
*/

double ikr_act, ikr_act_inf,tau_ikr_act;
double ikr;
double ikr_inact_inf,tau_ikr_inact,tau_ikr_inact2,ikr_inact,ikr_inact2;

/*
IKs variables
*/

double iks_act, iks_act_inf,tau_iks_act;
double iks;

/*
ICaL 1.2 and 1.3 parameters.
*/

double alpha_dl, beta_dl, tau_dl, alpha_fl, beta_fl, tau_fl;
double fl12, dl12, dl12_inf, fl12_inf,ical12;
double fl13, dl13, dl13_inf, fl13_inf,ical13;

double fca; 
double fca_inf ;
double taufca ;

/*
INa Nav1.1 (TTXS) and Nav1.5 (TTXR) variables,
*/

double ina_ttxr, ina_ttxs; // ina: Nav 1.1 and Nav 1.5
double m3_inf_ttxr, h_inf_ttxr;
double m3_inf_ttxs, h_inf_ttxs;
double m_inf_ttxr,m_inf_ttxs, hs,hsr;
double m_ttxr,h_ttxr,j_ttxr;
double m_ttxs,h_ttxs,j_ttxs;
double tau_m,tau_h,tau_j,tau_mr,tau_hr,tau_jr;
double delta_ttxr_m,delta_ttxr_h,delta_ttxr_j;
double fna;

/*
If variables
*/

double ih, ihk, ihna, ih_1, ih_2, ih_4, y_1_2, y_4, tau_y_1_2, tau_y_4, y_inf, y_1_2_inf, y_4_inf,ih_1_k, ih_2_k, ih_4_k, ih_1_na, ih_2_na, ih_4_na;

/*
INaK variables
*/

double inak;

/*
Inaca variables.
*/
double inaca;

double di,doo,k43,k12,k14,k41,k34,k21,k23,k32,x1,x2,x3,x4;

/*
Isus variables
*/

double r_inf, r, tau_r, isus;

/*
Ito variables
*/

double q,q_inf,tau_q;
double ito;

/*
Ionic homeostasis
*/

double ca_flux;

/*
Calcium diffusion
*/

double Jcadif;

double Jrel, Jup, Jtr,carel,caup;

double variation_multiply = 4.0;

double dFtc,Ftc;
double dFtmc,Ftmc,Ftmm;
double dFtmm,Ftmm;
double dFcms,Fcms;
double dFcmi,Fcmi;
double dFcq,Fcq;

int param_counter = 0, number = 0;
double min_potential[number_of_apds];
double tmin_potential[number_of_apds];
double max_potential[number_of_apds];
double dvdtmax[number_of_apds];
double vdvdtmax[number_of_apds];
double apd_start[number_of_apds], apd_end[number_of_apds],cai_peak[number_of_apds],cai_min[number_of_apds];
double ddr[number_of_apds];
double top[number_of_apds];
double top_slope[number_of_apds];
double apd50[number_of_apds];
double apd90[number_of_apds];
double cycle_length[number_of_apds];
double casub_peak[number_of_apds],casub_min[number_of_apds];

double dummy = 1.0;
double base_cycle = 0.0;
double dvdt, dvdtold;
double nai_tot, ki_tot;
double nai_integral, ki_integral;
double ki_integral1, ki_integral2;
double nai_min, nai_max, ki_min, ki_max;

double nai=8.0;         // mM 
double ki=140;
double FRT, RTF;

/*
RyR markov chain
*/

double resting, open, inactivated, resting_inactivated;
double kcasr,kosrca,kisrca;

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


/*
Perturbation parameters.
*/

double Pup;
double ks;
double gna_ttxs;
double gna_ttxr;
double gcat;
double gkr;
double gto;
double gcal12;
double gcal13;
double gh;
double vhalf_gh;
double kNaCa;

int  gcal_gcat = 0;

int broken_counter = 0;

int main(int count, char *argv[])
{  
FILE *outputcurrents,*parameters, *random;

char *str;

int filecounter = 0, start_output = 0, start_output1 = 0, filecounter1 = 1;
int filecountermin, filecountermax,filecounter1min, filecounter1max;

double total_time = 0.0;

int temp_counter = 0;
int one_or_two_pars = 0;

int param_basal;
int param_gna11;
int param_gna15;
int param_gcal12;
int param_gcal13;
int param_gcat;
int param_gf;
int param_vhalf_if;
int param_gkr;
int param_gto;
int param_knaca;
int param_pup;
int param_ks;

if(screen_output==1)
printf("%d \n",count);

/*
Initialising the random number generator is important.
*/

random = fopen("/dev/urandom","rb");
fread(&somerandomnumber, sizeof(unsigned int), 1, random);
fclose(random);

init_genrand(somerandomnumber);


/*
initialise all these parameters.
*/

 param_basal = 0;

 param_gna11 = 0;
 param_gna15 = 0;
 param_gcal12 = 0;
 param_gcal13 = 0;
 param_gcat = 0;
 param_gf = 0;
 param_vhalf_if = 0;
 param_gkr = 0;
 param_gto = 0;
 param_knaca = 0;
 param_pup = 0;
 param_ks = 0;

for(temp_counter = 1; temp_counter < count;temp_counter++){

if(atoi(argv[temp_counter])<0||atoi(argv[temp_counter])>1){ // this may still take floats, so just be careful.
 printf("bad input, exiting\n");
 exit(0);
}
}

if(count!=13){ printf("it takes integers as input\n"); exit(0);}

if(count==13){

 param_gna11 = atoi(argv[1]);
 param_gna15 = atoi(argv[2]);
 param_gcal12 = atoi(argv[3]);
 param_gcal13 = atoi(argv[4]);
 param_gcat = atoi(argv[5]);
 param_gf = atoi(argv[6]);
 param_vhalf_if = atoi(argv[7]);
 param_gkr = atoi(argv[8]);
 param_gto = atoi(argv[9]);
 param_knaca = atoi(argv[10]);
 param_pup = atoi(argv[11]);
 param_ks = atoi(argv[12]);

}

if(screen_output==1)
printf("%s %d %d %d %d %d %d %d %d %d %d %d %d\n",argv[0],param_gna11,
 param_gna15,
 param_gcal12,
 param_gcal13,
 param_gcat,
 param_gf,
 param_vhalf_if,
 param_gkr,
 param_gto,
 param_knaca,
 param_pup,
 param_ks);

// set values of the runs.

if(
 param_gna11==0&&
 param_gna15==0&&
 param_gcal12==0&&
 param_gcal13==0&&
 param_gcat==0&&
 param_gf==0&&
 param_vhalf_if==0&&
 param_gkr==0&&
 param_gto==0&&
 param_knaca==0&&
 param_pup==0&&
 param_ks==0){

filecounter1min = 0;
filecounter1max = 1;

filecountermin = 0;
filecountermax = 16;
}
else
{
filecounter1min = 0;
filecounter1max = 5500;

filecountermin = 0;
filecountermax = 0;
}

one_or_two_pars = 0;

one_or_two_pars=param_gna11+param_gna15+param_gcal12+param_gcal13+param_gcat+param_gf+param_vhalf_if+param_gkr+param_gto+param_knaca+param_pup+param_ks;

if(one_or_two_pars==2)
variation_multiply = 2.0;
else
if(one_or_two_pars==1)
variation_multiply = 4.0;
else
if(one_or_two_pars==0)
variation_multiply = 1.0;
else{
printf("this is bad\n");
exit(1);
}

FRT = F/(R*T);

broken_counter = 0;

for(filecounter1 = filecounter1min; filecounter1 < filecounter1max;filecounter1++){
for(filecounter = filecountermin; filecounter <=filecountermax; filecounter++)
{

param_counter = 0;
current_trace_counter = 0;

nai_integral = 0.0;
ki_integral = 0.0;
ki_integral1 = 0.0;
ki_integral2 = 0.0;

if(if_ap_profiles==1){
 str = malloc(32*sizeof(char));
 sprintf(str,"current%d_%d.dat",filecounter,filecounter1);
 outputcurrents = fopen(str,"w");
 free(str);
 }

/*
Initial conditions here.
*/

/* non-zero initial conditions required because 0 is a stable e.p. of the oscillator. */

resting = 0.74;
open = 0.0000034;
inactivated = 0.0000011;
resting_inactivated = 0.25;

nai_tot = 0.0;
ki_tot = 0.0;

v = -64.6507859122;
dst = 0.6185744228;
fst = 0.4571807567;
dt = 0.0015910667;
ft = 0.4307931823;
ikr_act = 0.4033909565;
ikr_inact = 0.9254741869;
ikr_inact2 = 0.1875749806;
iks_act = 0.0127072502;
fl12 = 0.9969251137;
dl12 = 0.0000044417;
fl13 = 0.9812942363;
dl13 = 0.0001993204;
r = 0.0045588036;
m_ttxr = 0.3991766512;
h_ttxr = 0.2765302400;
j_ttxr = 0.0252680441;
m_ttxs = 0.1069854109;
h_ttxs = 0.4548015654;
j_ttxs = 0.0272299711;
y_1_2 = 0.0281230987;
y_4 = 0.0137659036;
carel = 0.1038639395;
caup = 1.4120022486;
casub = 0.0000534745;
Ftc = 0.0065386464;
Ftmc = 0.1771999988;
Ftmm = 0.7268211394;
Fcms = 0.0231295949;
Fcmi = 0.0142601968;
Fcq = 0.1067110067;
cai = 0.0000326893;
q = 0.6140712332;
fca = 0.7679505412;
nai = 8.0839673080;
ki = 139.9228278979;
resting = 0.7718214188;
open = 0.0000000737;
inactivated = 0.0000000207;
resting_inactivated = 0.2173201024;

/*
Basal values of important parameters.
*/

gna_ttxs = 0.1*5.925e-05;
gna_ttxr = 0.1*5.925e-05;
gcal12 = 0.0010*4.0*1.5;
gcal13 = 0.0030*4.0*1.5;
gcat = 0.75*0.01862; 
gh = 0.0057; 
vhalf_gh = 106.8;
gkr = 0.8*0.002955;
gto = 0.00492;
 kNaCa = 5.5;
Pup = 0.04;
ks = 1300000.0;

/* genrand_real1(); random numbers between 0 and 1 inclusive. */
// 1 
if(param_gna11==1)
gna_ttxs = genrand_real1()*variation_multiply*0.1*5.925e-05; /* 0.05 x e-5 to 6e-5 */
// 2
if(param_gna15==1)
gna_ttxr = genrand_real1()*variation_multiply*0.1*5.925e-05;
// 3
if(param_gcal12==1)
gcal12 = genrand_real1()*variation_multiply*0.0010*4.0*1.5;
// 4
if(param_gcal13==1)
gcal13 = genrand_real1()*variation_multiply*0.0030*4.0*1.5; /* 1 to 12 multiply */
// 5
if(param_gcat==1)
gcat = genrand_real1()*variation_multiply*0.75*0.01862; /* 0.005 to 0.2 */
// 6
if(param_gf==1)
gh = genrand_real1()*variation_multiply*0.0057; /* 0 to 3 times */
// 7
if(param_vhalf_if==1)
vhalf_gh = genrand_real1()*variation_multiply*106.8; // between -160 and -50 mV
// 8
if(param_gkr==1)
gkr = genrand_real1()*variation_multiply*0.8*0.002955;      /* 0.002955 to 5*0.002955 */
// 9
if(param_gto==1)
gto = genrand_real1()*variation_multiply*0.00492;     /* 0.0005 to 0.01 */
// 10
if(param_knaca==1)
kNaCa = genrand_real1()*variation_multiply*5.5; /* 4 to 12 */
// 11
if(param_pup==1)
Pup = genrand_real1()*variation_multiply*0.04; /* 0.01 to 0.3 */
// 12
if(param_ks==1)
ks = genrand_real1()*variation_multiply*1300000.0; /* 0 to 1800000 */

gcal_gcat = 0;
if(gcat<=gcal13){
gcal_gcat = 1; // this doesnt mean much, gcat is actually gcatmax*act*inact (t) or ever time averaged/clamped.
}

for(param_counter=0;param_counter<number_of_apds;param_counter++){
     min_potential[param_counter] = 10000.0; // impossibly large.
     max_potential[param_counter] = -10000.0;
     dvdtmax[param_counter] = -10000.0;
     cai_peak[param_counter] = -1000000.0;
     cai_min[param_counter] = 10000.0;
     ddr[param_counter] = -10000.0;
     top[param_counter] = 10000.0;
     top_slope[param_counter] = -10000.0;
     apd50[param_counter] = -10000.0;
     apd90[param_counter] = -10000.0;
     cycle_length[param_counter] = -10000.0;
     casub_peak[param_counter] = -1000000.0;
     casub_min[param_counter] = 10000.0;
     
     if(param_counter<=filecountermax)     
	base_cycle = -1.0;
	
	solution_broke = 0;
}

/*
Reversal potentials
*/

ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);

param_counter = 0;

dvdt = 1000.0;
dvdtold = 500.0;

start_output = 0;
start_output1 = 0;

total_time = 12000.0;

for(sstime=0.0;sstime<=total_time;sstime = sstime + ddt){         // time loop starts here

if(v!=v){
printf("solution broke %d!\n",filecounter);
break;
}

/*Ist********************************************************************/
/* Shinagawa  */

qa = 1.0/(1.0 + exp(-(v+67.0)/5.0));

alphaqa = 1.0/(0.15*exp(-(v)/11.0)+0.2*exp(-(v)/700.0));
betaqa  =  1.0/(16.0*exp((v)/8.0)+15.0*exp((v)/50.0));
tauqa = 1.0/(alphaqa + betaqa);

alphaqi = 0.15*1.0/(3100.0*exp((v+10.0)/13.0)+700.3*exp((v+10.0)/70.0));
betaqi =  0.15*1.0/(95.7*exp(-(v+10.0)/10.0) + 50.0*exp(-(v+10.0)/700.0)) + 0.000229/(1+exp(-(v+10.0)/5.0));
qi = alphaqi/(alphaqi + betaqi);
tauqi = 1.0/(alphaqi + betaqi);

dst = dst + ddt*((qa-dst)/tauqa);
fst = fst + ddt*((qi-fst)/tauqi);

if((filecounter==6)&&(whichParamForBif==0))
ist = 0.0*gst*dst*fst*(v - eist);
else
ist = gst*dst*fst*(v - eist);

/* Ib ************************************************************************/
  ibna = gbna*(v - ena);
  ibca = gbca*(v - eca);
  ibk  =  gbk*(v - ek);
  ib = (ibna + ibca + ibk);
  
/*IK1**********************************************************************/
xk1inf = 1.0/(1.0 + exp(0.070727*(v - ek)));
ik1 = gk1*xk1inf*(ko/(ko + 0.228880))*(v - ek);

/**ICaT Cav3.1**************************************************************/
/* ICaT */

tau_dt = 1.0/(1.068*exp((v + 26.3)/30.0) + 1.068*exp(-(v + 26.3)/30.0)); // Kurata
dt_inf = 1.0/(1.0+exp(-(v + 26.0)/6.0));
dt = dt + ddt*((dt_inf - dt)/tau_dt);

tau_ft = 1.0/(0.0153*exp(-(v+61.7)/83.3)+0.015*exp((v+61.7)/15.38)); // Kurata/Maltsev.
ft_inf = 1.0/(1.0+exp((v + 61.7)/5.6)); // rabbit SS, since the Mangoni data is 10 mV low, and very steep.
ft = ft + ddt*((ft_inf - ft)/tau_ft);

if((filecounter==7)&&(whichParamForBif==0))
icat = 0.0*gcat*ft*dt*(v - ecat);           // nA
else
icat = gcat*ft*dt*(v - ecat);           // nA

/*Ikr*************************************************************************/

if((filecounter==13)&&(whichParamForBif==0))
ikr_act_inf = 1.0/(1.0 + exp(-(v+21.173694+5.0)/9.757086)); // ISO shifts it to the left.
else
ikr_act_inf = 1.0/(1.0 + exp(-(v+21.173694)/9.757086)); // sacred
tau_ikr_act = 0.699821/(0.003596*exp((v)/15.339290) + 0.000177*exp(-(v)/25.868423)); // this fits the Q10 variation.
ikr_act = ikr_act + ddt*(ikr_act_inf-ikr_act)/tau_ikr_act;
     
/* SS inactivation */

ikr_inact_inf = 1.0/(1.0 + exp((v+20.758474-4.0)/(19.0))); 

// ikr_inact_inf = 1.0/(1.0 + exp((v+28.6)/(17.1)));  // this was choosen in Jan 2011.
// ikr_inact_inf = 1.0/(1.0 + exp((v+15.5)/5.3)); // mangoni
// ikr_inact_inf = 1.0/(1.0 + exp((v+20.7)/5.3));  // sk mod

tau_ikr_inact = 0.2+0.9*1.0/(0.1*exp(v/54.645)+0.656*exp(v/106.157));

//tau_ikr_inact = 0.84655/(0.0372*exp(v/15.9)+0.00096*exp(v/22.5)); // the fast one
//tau_ikr_inact2 = 0.84655/(0.0042*exp(v/17.0)+0.00015*exp(v/21.6)); // the slow one

ikr_inact = ikr_inact + ddt*(ikr_inact_inf - ikr_inact)/tau_ikr_inact;
// ikr_inact2 = ikr_inact2 + ddt*(ikr_inact_inf - ikr_inact2)/tau_ikr_inact2; // no evidence of biphasic activation, or inactivation.

if((filecounter==8)&&(whichParamForBif==0))
ikr = 0.0*gkr*ikr_act*ikr_inact*(v - ek);
else
if((filecounter==13)&&(whichParamForBif==0))
ikr = 1.12*gkr*ikr_act*ikr_inact*(v - ek);
else
if((filecounter==15)&&(whichParamForBif==0))
ikr = 0.52*gkr*ikr_act*ikr_inact*(v - ek);
else
ikr = gkr*ikr_act*ikr_inact*(v - ek);

/**IKs********************************************************************/

iks_act_inf = 1.0/(1.0 + exp(-(v-20.876040)/11.852723));
tau_iks_act =  1000.0/(13.097938/(1.0 + exp(-(v-48.910584)/10.630272)) + exp(-(v)/35.316539)); // Zhang model. The other models are similar.
iks_act = iks_act + ddt*(iks_act_inf - iks_act)/tau_iks_act; // Ding/Maatsura
iks = gks*iks_act*iks_act*(v - eks);

/*ICaL*******************************************************************/

if(fabs(v)<=0.001){
 alpha_dl  = -28.39*(v+35.0)/(exp(-(v+35.0)/2.5)-1.0)+408.173;
}
else
if(fabs(v+35.0)<=0.001){
 alpha_dl  = 70.975-84.9*v/(exp(-0.208*v)-1.0);
}
else
if(fabs(v)>0.001&&fabs(v+35.0)>0.001){
 alpha_dl  = -28.39*(v+35.0)/(exp(-(v+35.0)/2.5)-1.0)-84.9*v/(exp(-0.208*v)-1.0);
}

if(fabs(v-5.0)<=0.001)
 beta_dl   = 28.575;
else
if(fabs(v-5.0)>0.001)
 beta_dl   = 11.43*(v-5.0)/(exp(0.4*(v-5.0))-1.0);

tau_dl  = 2000.0/(alpha_dl +beta_dl);

/* Cav 1.3 */
/* Cav 1.3 SS values */

// dl13_inf = 1.0/(1+exp(-(v+25.0)/5.39)); // Mangoni 2006 who took it from Zhang
// fl13_inf = 1.0/(1+exp((v+45.66)/4.1));  // Mangoni 2006 who took it from Zhang

dl13_inf = 1.0/(1+exp(-(v+13.5)/6.0)); 
fl13_inf = 1.0/(1+exp((v+35.0)/7.3)); 

tau_fl = 7.4 + 45.77*exp(-0.5*(v+28.1)*(v+28.1)/(11*11));
// tau_fl = 44.3 + 257.1*exp(-(v+32.5)*(v+32.5)/(13.9*13.9)); // Maltsev TC. Its a bit off for our model, but it may fit if tried upon.

dl13 = dl13 + ddt*(dl13_inf - dl13)/tau_dl;
fl13 = fl13 + ddt*(fl13_inf - fl13)/tau_fl;

/* Cav 1.2 */

dl12_inf = 1.0/(1+exp(-(v+3.0)/5.0)); // Mangoni 2006: according to his email
fl12_inf = 1.0/(1+exp((v+36.0)/4.6)); // Mangoni 2006 35.9

dl12 = dl12 + ddt*(dl12_inf - dl12)/tau_dl;
fl12 = fl12 + ddt*(fl12_inf - fl12)/tau_fl;

/* fca */

fca_inf = kmfca/(kmfca+casub);
taufca = fca_inf/alpha_fca;
fca = fca + ddt*(fca_inf - fca)/taufca;

if((filecounter==9)&&(whichParamForBif==0))
 ical12 = 0.0*gcal12*fl12*dl12*fca*(v-ecal);
else
if((filecounter==13)&&(whichParamForBif==0))
 ical12 = 1.45*gcal12*fl12*dl12*fca*(v-ecal);
else
 ical12 = gcal12*fl12*dl12*fca*(v-ecal);
    
if((filecounter==10)&&(whichParamForBif==0))
 ical13 = 0.0*gcal13*fl13*dl13*fca*(v-ecal);
else
if((filecounter==13)&&(whichParamForBif==0))
 ical13 = 1.45*gcal13*fl13*dl13*fca*(v-ecal);
else
if((filecounter==14)&&(whichParamForBif==0))
 ical13 = 0.56*gcal13*fl13*dl13*fca*(v-ecal);
else
 ical13 = gcal13*fl13*dl13*fca*(v-ecal);
    
/**INa**************************************************************************/

fna = (9.52e-02*exp(-6.3e-2*(v+34.4))/(1+1.66*exp(-0.225*(v+63.7))))+8.69e-2; // 1-fna = proportion of fast inactivation, fna = contribution of slow

m3_inf_ttxr = 1.0/(1.0 + exp(-(v+45.213705)/7.219547));
h_inf_ttxr = 1.0/(1.0 + exp(-(v+62.578120 )/(-6.084036)));
// m3_inf_ttxs = 1.0/(1.0 + exp(-(v+36.097331)/6.258247));
m3_inf_ttxs = 1.0/(1.0 + exp(-(v+36.097331-5.0)/5.0)); // Mangoni takes -29.
h_inf_ttxs = 1.0/(1.0 + exp((v+56.0)/3.0));

m_inf_ttxr = pow(m3_inf_ttxr,0.333);
m_inf_ttxs = pow(m3_inf_ttxs,0.333);

tau_m = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_h = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_j = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);
		
m_ttxs = m_ttxs + ddt*(m_inf_ttxs - m_ttxs)/tau_m;
h_ttxs = h_ttxs + ddt*(h_inf_ttxs - h_ttxs)/tau_h;
j_ttxs = j_ttxs + ddt*(h_inf_ttxs - j_ttxs)/tau_j;
 
hs = (1.0-fna)*h_ttxs+fna*j_ttxs;

tau_mr = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_hr = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_jr = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);

m_ttxr = m_ttxr + ddt*(m_inf_ttxr - m_ttxr)/tau_mr;
h_ttxr = h_ttxr + ddt*(h_inf_ttxr - h_ttxr)/tau_hr;
j_ttxr = j_ttxr + ddt*(h_inf_ttxr - j_ttxr)/tau_jr;

hsr = (1.0-fna)*h_ttxr+fna*j_ttxr;

if((filecounter==11)&&(whichParamForBif==0)){
if(fabs(v)>0.005)
ina_ttxs=0.0*gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*(F*F/(R*T))*((exp((v-ena)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxs=0.0*gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*F*((exp((v-ena)*F/(R*T))-1.0));
}
else
{
if(fabs(v)>0.005)
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*(F*F/(R*T))*((exp((v-ena)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*F*((exp((v-ena)*F/(R*T))-1.0));
}

if((filecounter==12)&&(whichParamForBif==0)){
if(fabs(v)>0.005)
ina_ttxr = 0.0*gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*(F*F/(R*T))*((exp((v-enattxr)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxr = 0.0*gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*F*((exp((v-enattxr)*F/(R*T))-1.0));
}
else
{
if(fabs(v)>0.005)
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*(F*F/(R*T))*((exp((v-enattxr)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*F*((exp((v-enattxr)*F/(R*T))-1.0));
}

/**If**************************************************************************/

if((filecounter==13)&&(whichParamForBif==0)){
y_inf = 1.0/(1.0 + exp((v+92.0)/16.3));
}
else
y_inf = 1.0/(1.0 + exp((v+vhalf_gh)/16.3));
tau_y_1_2 = 1.5049/(exp(-(v+590.3)*0.01094)+ exp((v-85.1)/17.2));
y_1_2 = y_1_2 + ddt*(y_inf - y_1_2)/tau_y_1_2;

if((filecounter==5)&&(whichParamForBif==0)){
ihk = 0.0;
ihna = 0.0;
 }
else
{
ihk = 0.6167*gh*y_1_2*(v - ek);
ihna = 0.3833*gh*y_1_2*(v - ena);
}

ih = (ihk + ihna);

/*Ito*************************************************************************/

q_inf = 1.0/(1.0+exp((v+49.0)/13.0));
tau_q = (6.06 + 39.102/(0.57*exp(-0.08*(v+44.0))+0.065*exp(0.1*(v+45.93))))/0.67; 
q = q + ddt*((q_inf-q)/tau_q);
r_inf = 1.0/(1.0+exp(-(v-19.3)/15.0));
tau_r = (2.75+14.40516/(1.037*exp(0.09*(v+30.61))+0.369*exp(-0.12*(v+23.84))))/0.303;
r = r + ddt*((r_inf-r)/tau_r);
ito = gto*q*r*(v-ek);
 
/*Isus***********************************************************************/

isus = gsus*r*(v-ek);

/*Inak***********************************************************************/

// Kurata/Maltsev INaK
 inak = inakmax*(pow(ko,1.2)/(pow(kmkp,1.2)+pow(ko,1.2)))*(pow(nai,1.3)/(pow(kmnap,1.3)+pow(nai,1.3)))/(1.0+exp(-(v-ena+120.0)/30.0));

/****iNaCa*******************************************************************/

// INaCa
       di=1+(casub/Kci)*(1+exp(-Qci*v*FRT)+nai/Kcni)+(nai/K1ni)*(1+(nai/K2ni)*(1+nai/K3ni));
       doo=1+(cao/Kco)*(1+exp(Qco*v*FRT))+(nao/K1no)*(1+(nao/K2no)*(1+nao/K3no));
       k43=nai/(K3ni+nai);
       k12=(casub/Kci)*exp(-Qci*v*FRT)/di;
       k14=(nai/K1ni)*(nai/K2ni)*(1+nai/K3ni)*exp(Qn*v*FRT/2.0)/di;
       k41=exp(-Qn*v*FRT/2.0);
       k34=nao/(K3no+nao);
       k21=(cao/Kco)*exp(Qco*v*FRT)/doo;
       k23=(nao/K1no)*(nao/K2no)*(1+nao/K3no)*exp(-Qn*v*FRT/2.0)/doo;
       k32=exp(Qn*v*FRT/2);
       x1=k34*k41*(k23+k21)+k21*k32*(k43+k41);
       x2=k43*k32*(k14+k12)+k41*k12*(k34+k32);
       x3=k43*k14*(k23+k21)+k12*k23*(k43+k41);
       x4=k34*k23*(k14+k12)+k21*k14*(k34+k32);
   
if((filecounter==3)&&(whichParamForBif==0))
 inaca = 0.0;
else
if((filecounter==16)&&(whichParamForBif==0))
 inaca = 0.3*kNaCa*(k21*x2-k12*x1)/(x1+x2+x3+x4);   
else
 inaca = kNaCa*(k21*x2-k12*x1)/(x1+x2+x3+x4);   

/*****************************************************************************/

/* Ion fluxes */
/* complete these equations */
ca_flux = (ical12+ical13+icat-2.0*inaca+ibca)/(2.0*F);

/*****************************************************************************/

 /* Ca2+ difusion */
 Jcadif = (casub - cai)/tdifca;

 /* Ca handling in SR (mM/ms) */
 kcasr = maxsr - (maxsr - minsr)/(1.0 + pow(eca50sr/carel,hsrr));
 kosrca = koca/kcasr;
 kisrca = kica*kcasr;

resting = resting + ddt*(kim*resting_inactivated - kisrca*casub*resting - kosrca*casub*casub*resting + kom*open);
open    = open    + ddt*(kosrca*casub*casub*resting - kom*open - kisrca*casub*open + kim*inactivated);
inactivated = inactivated + ddt*(kisrca*casub*open - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated);
resting_inactivated = resting_inactivated + ddt*(kom*inactivated - kosrca*casub*casub*resting_inactivated - kim*resting_inactivated + kisrca*casub*resting);

if((filecounter==4)&&(whichParamForBif==0))
Jrel = 0.0*ks*open*(carel - casub);
else
Jrel = ks*open*(carel - casub);

if((filecounter==2)&&(whichParamForBif==0))
 Jup  = 0.0*Pup/(1.0+Kup/cai);
else
if((filecounter==13)&&(whichParamForBif==0))
 Jup  = 4.0*Pup/(1.0+Kup/cai);    
else
 Jup  = Pup/(1.0+Kup/cai);
           
Jtr  = (caup - carel)/Ttr;

// Ca buffering flux // these are the derivatives of the F's
dFtc  = kfTC*cai*(1.0-Ftc)-kbTC*Ftc;
dFtmc = kfTMC*cai*(1.0-Ftmc-Ftmm)-kbTMC*Ftmc;
dFtmm = kfTMM*Mgi*(1.0-Ftmc-Ftmm)-kbTMM*Ftmm;
dFcms = kfCM*casub*(1.0-Fcms)-kbCM*Fcms;
dFcmi = kfCM*cai*(1.0-Fcmi)-kbCM*Fcmi;
dFcq  = kfCQ*carel*(1.0-Fcq)-kbCQ*Fcq;

Ftc = Ftc + ddt*dFtc;                   //	    % dfTnCa/dt
Ftmc = Ftmc + ddt*dFtmc;                  //	% dfTnMgCa/dt
Ftmm = Ftmm + ddt*dFtmm;                  //	% dfTnMgMg/dt
Fcms = Fcms + ddt*dFcms;                 //	    % dfCMs/dt
Fcmi = Fcmi + ddt*dFcmi;                 //	    % dfCMi/dt
Fcq = Fcq + ddt*dFcq;                    //	    % dfCQ/dt

if((filecounter==1)&&(whichParamForBif==0)){
 casub = 0.000000116; // casub + 0.0*ddt*((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);
 cai = 0.000000116; // cai + 0.0*ddt*((Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmi + ConcTC*dFtc + ConcTMC*dFtmc));
}
else
{
 casub = casub + ddt*((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);
 cai = cai + ddt*((Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmi + ConcTC*dFtc + ConcTMC*dFtmc));
}

carel = carel + ddt*(Jtr - Jrel - ConcCQ*dFcq);
caup = caup + ddt*(Jup-Jtr*vrel/vup);

/*************************************************************************************************/

dvdtold = dvdt;
total_current = ih+ina_ttxr+ina_ttxs+ical12+ical13+iks+ikr+ik1+ist+ib+icat+inak+isus+inaca+ito;
dvdt = - total_current/capacitance;

vnew = v  + ddt*dvdt;
 
ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);
 
nai_tot = ihna+ina_ttxr+ina_ttxs+3.0*inak+3.0*inaca+ist+ibna;
ki_tot = ihk+iks+ikr+ik1+ibk-2.0*inak+isus+ito;


if((filecounter==1)&&(whichParamForBif==0)){
 nai = nai;
 ki = ki;
}
else
if((filecounter==2)&&(whichParamForBif==0)){
 nai = nai;
 ki = ki;
}
else
if((filecounter==4)&&(whichParamForBif==0)){
 nai = nai;
 ki = ki;
}
else
{
nai = nai - ddt*(nai_tot)/(F*vi);
ki = ki - ddt*(ki_tot)/(F*vi);
}
 
nai_integral = nai_integral + ddt*(nai_tot)/(F*vi);
ki_integral = ki_integral + ddt*(ki_tot)/(F*vi);
ki_integral1 = ki_integral1 + ddt*(ihk+iks+ikr+ik1+ibk+isus+ito);
ki_integral2 = ki_integral2 + ddt*(-2.0*inak/inakmax);

/*
Measurements here.
*/

if(sstime>(total_time*0.1)){
if(cai_min[param_counter]>cai) cai_min[param_counter] = cai;
if(cai_peak[param_counter]<cai) cai_peak[param_counter] = cai;
if(casub_min[param_counter]>casub) casub_min[param_counter] = casub;
if(casub_peak[param_counter]<casub) casub_peak[param_counter] = casub;

 if(dvdt>=0.0&&dvdtold<0.0){
  min_potential[param_counter] = v;
  nai_min = nai;
  ki_min = ki;
  tmin_potential[param_counter] = sstime;
  start_output = 1;
  start_output1 = 1;
 }

if(dvdt>dvdtmax[param_counter]&&start_output>0){
 dvdtmax[param_counter] = dvdt;
 apd_start[param_counter] = sstime;
 vdvdtmax[param_counter] = v;	
}

if(dvdtold>0.0&&dvdt<=0.0){
 max_potential[param_counter] = v;
 nai_max = nai;
 ki_max = ki;
 top_slope[param_counter] = (max_potential[param_counter]-min_potential[param_counter])/(sstime - tmin_potential[param_counter]);
}

if((param_counter>0)&&(dvdtold<=top_slope[param_counter-1])&&(dvdt>top_slope[param_counter-1])){
 top[param_counter] = v;
 ddr[param_counter] = (v - min_potential[param_counter])/(sstime - tmin_potential[param_counter]);
}

if(vnew<=0.5*min_potential[param_counter]&&v>0.5*min_potential[param_counter]){
 if(apd_start[param_counter]>0.0)
 apd50[param_counter] = sstime - apd_start[param_counter];
}

if(vnew<=0.9*min_potential[param_counter]&&v>0.9*min_potential[param_counter]){
 if(apd_start[param_counter]>0.0){
   apd_end[param_counter] = sstime;
   apd90[param_counter] = sstime - apd_start[param_counter];
   cycle_length[param_counter] = apd_start[param_counter]-apd_start[param_counter-1];

if(screen_output==1){
     printf("%d\t",filecounter);
     printf("%10.5f\t",min_potential[param_counter]);
     printf("%10.5f\t",max_potential[param_counter]);
     printf("%10.5f\t",dvdtmax[param_counter]);
     printf("%10.5f\t",apd50[param_counter]);
     printf("%10.5f\t",apd90[param_counter]);
     printf("%10.5f\t",cycle_length[param_counter]);
     printf("%10.5f\t",ddr[param_counter]);
     printf("%10.5f\t",top[param_counter]);
     printf("%10.5f\t",casub_min[param_counter]);
     printf("%10.5f\t",casub_peak[param_counter]);
     printf("%10.5f\t",ki);
     printf("%10.5f\t",nai);
     printf("%20.10f\t",nai_integral);
     printf("%20.10f\t",ki_integral);
     printf("\n");
     }
     
     param_counter++;
     
     nai_integral = 0.0;
     ki_integral = 0.0;	
     ki_integral1 = 0.0;
     ki_integral2 = 0.0;
     }
   }
} // end of measurement

v = vnew;

/************************************************************************/

if(if_ap_profiles==1)
if(current_trace_counter%800==0&&sstime>(total_time-2000)){
// if(current_trace_counter%800==0&&sstime>0.05*total_time){
// 1 2
fprintf(outputcurrents,"%10.10f\t%10.10f\t",sstime-(total_time-2000),v);
// 3
fprintf(outputcurrents,"%10.10f\t",ina_ttxs/capacitance);
// 4
fprintf(outputcurrents,"%10.10f\t",ical12/capacitance);
// 5
fprintf(outputcurrents,"%10.10f\t",ical13/capacitance);
// 6
fprintf(outputcurrents,"%10.10f\t",iks/capacitance);
// 7
fprintf(outputcurrents,"%10.10f\t",ikr/capacitance);
// 8
fprintf(outputcurrents,"%10.10f\t",ih/capacitance);
// 9
fprintf(outputcurrents,"%10.10f\t",ik1/capacitance);
// 10
fprintf(outputcurrents,"%10.10f\t",ist/capacitance);
// 11
fprintf(outputcurrents,"%10.10f\t",ito/capacitance);
// 12
fprintf(outputcurrents,"%10.10f\t",isus/capacitance);
// 13
fprintf(outputcurrents,"%10.10f\t",inaca/capacitance);
// 14
fprintf(outputcurrents,"%10.10f\t",inak/capacitance);
// 15
fprintf(outputcurrents,"%10.10f\t",ib/capacitance);
// 16
fprintf(outputcurrents,"%10.10f\t",cai);
// 17
fprintf(outputcurrents,"%10.10f\t",total_current/capacitance);
// 18
fprintf(outputcurrents,"%10.10f\t",dvdt);
// 19
fprintf(outputcurrents,"%10.10f\t",icat/capacitance);
// 20
fprintf(outputcurrents,"%10.10f\t",ki);
// 21
fprintf(outputcurrents,"%10.10f\t",nai);
// 22
fprintf(outputcurrents,"%10.10f\t",ina_ttxr/capacitance);
// 23
fprintf(outputcurrents,"%10.10f\t",casub);
// 24
fprintf(outputcurrents,"%10.10f\t",caup);
// 25
fprintf(outputcurrents,"%10.10f\t",carel);
// 26
fprintf(outputcurrents,"%10.10f\t",Jrel);
// 27
fprintf(outputcurrents,"%10.10f\t",Jup);
// 28
fprintf(outputcurrents,"%10.10f\t",Jtr);
// 29
fprintf(outputcurrents,"%10.10f\t",Ftc);
// 30
fprintf(outputcurrents,"%10.10f\t",Ftmc);
// 31
fprintf(outputcurrents,"%10.10f\t",Ftmm);
// 32
fprintf(outputcurrents,"%10.10f\t",Fcms);
// 33
fprintf(outputcurrents,"%10.10f\t",Fcmi);
// 34
fprintf(outputcurrents,"%10.10f\t",Fcq);
// 35
fprintf(outputcurrents,"%10.10f\t",ical12/capacitance+ical13/capacitance);
// 36
fprintf(outputcurrents,"%10.10f\t",open);
// 37
fprintf(outputcurrents,"%10.10f\t",resting);
// 38
fprintf(outputcurrents,"%10.10f\t",inactivated);
// 39
fprintf(outputcurrents,"%10.10f\t",resting_inactivated);

fprintf(outputcurrents,"\n");
//current_trace_counter = 0;
}

// new initial conditions to reduce run times.
if(filecounter==0&&start_output1==1){
parameters = fopen("initial_conditions.dat","w");
fprintf(parameters,"v = %10.10f;\n",v); // v = -56.391127 ;
fprintf(parameters,"dst = %10.10f;\n",dst); // dst = 0.705069;
fprintf(parameters,"fst = %10.10f;\n",fst); // fst = 0.178856;
fprintf(parameters,"dt = %10.10f;\n",dt); // dt = 0.645266;
fprintf(parameters,"ft = %10.10f;\n",ft); // ft = 0.002712;
fprintf(parameters,"ikr_act = %10.10f;\n",ikr_act); // ikr_act = 0.556028;
fprintf(parameters,"ikr_inact = %10.10f;\n",ikr_inact); // ikr_inact = 0.951085;
fprintf(parameters,"ikr_inact2 = %10.10f;\n",ikr_inact2); // ikr_inact = 0.951085;
fprintf(parameters,"iks_act = %10.10f;\n",iks_act); // iks_act = 0.104055;
fprintf(parameters,"fl12 = %10.10f;\n",fl12); // fl12 = 0.814427;
fprintf(parameters,"dl12 = %10.10f;\n",dl12); // dl12 = 0.000635;
fprintf(parameters,"fl13 = %10.10f;\n",fl13); // fl13 = 0.535972;
fprintf(parameters,"dl13 = %10.10f;\n",dl13); // dl13 = 0.057282;
fprintf(parameters,"r = %10.10f;\n",r); // r = 0.044534;
fprintf(parameters,"m_ttxr = %10.10f;\n",m_ttxr); // m_ttxr = 0.552195;
fprintf(parameters,"h_ttxr = %10.10f;\n",h_ttxr); // h_ttxr = 0.064120;
fprintf(parameters,"j_ttxr = %10.10f;\n",j_ttxr); // j_ttxr = 0.115122;
fprintf(parameters,"m_ttxs = %10.10f;\n",m_ttxs); // m_ttxs = 0.501233;
fprintf(parameters,"h_ttxs = %10.10f;\n",h_ttxs); // h_ttxs = 0.110766;
fprintf(parameters,"j_ttxs = %10.10f;\n",j_ttxs); // j_ttxs = 0.181200;
fprintf(parameters,"y_1_2 = %10.10f;\n",y_1_2); // y_1_2 = 0.021433;
fprintf(parameters,"y_4 = %10.10f;\n",y_4); // y_4 = 0.008980;
fprintf(parameters,"carel = %10.10f;\n",carel); // carel = 8.943277;
fprintf(parameters,"caup = %10.10f;\n",caup); // caup = 10.159016;
fprintf(parameters,"casub = %10.10f;\n",casub); // casub = 0.000239;
fprintf(parameters,"Ftc = %10.10f;\n",Ftc); // Ftc = 0.267413;
fprintf(parameters,"Ftmc = %10.10f;\n",Ftmc); // Ftmc = 0.872259;
fprintf(parameters,"Ftmm = %10.10f;\n",Ftmm); // Ftmm = 0.125869;
fprintf(parameters,"Fcms = %10.10f;\n",Fcms); // Fcms = 0.384077;
fprintf(parameters,"Fcmi = %10.10f;\n",Fcmi); // Fcmi = 0.544780;
fprintf(parameters,"Fcq = %10.10f;\n",Fcq); // Fcq = 1.139253;
fprintf(parameters,"cai = %10.10f;\n",cai); // cai = 0.000097;
fprintf(parameters,"q = %10.10f;\n",q); // q = 0.255269;
fprintf(parameters,"fca = %10.10f;\n",fca); // fca = 0.550188;
fprintf(parameters,"nai = %10.10f;\n",nai); // nai=8.0;         
fprintf(parameters,"ki = %10.10f;\n",ki); // ki=140;       
fprintf(parameters,"resting = %10.10f;\n",resting); // ki=140;       
fprintf(parameters,"open = %10.10f;\n",open); // ki=140;       
fprintf(parameters,"inactivated = %10.10f;\n",inactivated); // ki=140;       
fprintf(parameters,"resting_inactivated = %10.10f;\n",resting_inactivated); // ki=140; 
fclose(parameters);

start_output1=0;
}

current_trace_counter++;

} // end of time loop

if(if_ap_profiles==1)
fclose(outputcurrents);

str = malloc(32*sizeof(char));
// sprintf(str,"APparams%d.dat",filecounter1); // making sure only sane data is outputted. if insane data, the file is overwritten.
sprintf(str,"APparams.dat");
number = param_counter;
parameters = fopen(str,"a+");
 free(str);

if((one_or_two_pars==2)&&(number>=3)){
for (param_counter=number-3; param_counter< number; param_counter++)
fprintf( parameters,"%d\t%d\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%d\n",
filecounter,filecounter1,
gna_ttxs,gna_ttxr,
gcal12,gcal13,
gcat,gh,
vhalf_gh,gkr,
gto,kNaCa,
Pup,ks,
cycle_length[param_counter],-min_potential[param_counter]+max_potential[param_counter],
min_potential[param_counter],max_potential[param_counter],
dvdtmax[param_counter],apd50[param_counter],
apd90[param_counter],cai_min[param_counter],
cai_peak[param_counter],ddr[param_counter],
top[param_counter],gcal_gcat);
}
else
if((one_or_two_pars==1)&&(number>=4)){
for (param_counter=number-4; param_counter< number; param_counter++)
fprintf( parameters,"%d\t%d\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%d\n",
filecounter,filecounter1,
gna_ttxs,gna_ttxr,
gcal12,gcal13,
gcat,gh,
vhalf_gh,gkr,
gto,kNaCa,
Pup,ks,
cycle_length[param_counter],-min_potential[param_counter]+max_potential[param_counter],
min_potential[param_counter],max_potential[param_counter],
dvdtmax[param_counter],apd50[param_counter],
apd90[param_counter],cai_min[param_counter],
cai_peak[param_counter],ddr[param_counter],
top[param_counter],gcal_gcat);
}
else
if((one_or_two_pars==0)&&(number>=4)){
for (param_counter=number-4; param_counter< number; param_counter++)
fprintf( parameters,"%d\t%d\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%d\n",
filecounter,filecounter1,
gna_ttxs,gna_ttxr,
gcal12,gcal13,
gcat,gh,
vhalf_gh,gkr,
gto,kNaCa,
Pup,ks,
cycle_length[param_counter],-min_potential[param_counter]+max_potential[param_counter],
min_potential[param_counter],max_potential[param_counter],
dvdtmax[param_counter],apd50[param_counter],
apd90[param_counter],cai_min[param_counter],
cai_peak[param_counter],ddr[param_counter],
top[param_counter],gcal_gcat);
}
else
{
fprintf( parameters,"%d\t%d\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%d\n",
filecounter,filecounter1,gna_ttxs,gna_ttxr,gcal12,gcal13,gcat,gh,vhalf_gh,gkr,gto,kNaCa,Pup,ks,-100.0,-100.0,100.0,-100.0,-100.0,-100.0,-100.0,-100.0,-100.0,-100.0,-100.0,gcal_gcat);
}		   
fclose(parameters);

} // end of filecounter loop.
} // end of filecounter1 loop. 

return 0;
}
