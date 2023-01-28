/*
* f routine. Compute function f(t,y). 
*/

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{

 realtype Y[NEQS];
 realtype dY[NEQS];
 int i;
 realtype IV;
 UserData data;

  for(i=0;i<NEQS;i++)
    Y[i] = Ith(y,i+1);

  data = (UserData) user_data;

/*data->p[0] =*/ realtype    R ; // data->p[0];
/*data->p[1] =*/ realtype    T ; // data->p[1];
/*data->p[2] =*/ realtype    F ; // data->p[2];
/*data->p[3] =*/ realtype    capacitance ; // data->p[3];
/*data->p[4] =*/ realtype    vcell ; // data->p[4];
/*data->p[5] =*/ realtype    l_cell ; // data->p[5];
/*data->p[6] =*/ realtype    r_cell ; // data->p[6];
/*data->p[7] =*/ realtype    vrel ; // data->p[7];
/*data->p[8] =*/ realtype    vsub ; // data->p[8];
/*data->p[9] =*/ realtype    vup  ; // data->p[9];
/*data->p[10] =*/ realtype    vi ; // data->p[10];
/*data->p[11] =*/ realtype    Mgi ; // data->p[11];
/*data->p[12] =*/ realtype    nao ; // data->p[12];
/*data->p[13] =*/ realtype    cao ; // data->p[13];
/*data->p[14] =*/ realtype    ko ; // data->p[14];
/*data->p[15] =*/ realtype    gst ; // data->p[15];
/*data->p[16] =*/ realtype    eist ; // data->p[16];
/*data->p[17] =*/ realtype    gbna ; // data->p[17];
/*data->p[18] =*/ realtype    gbca ; // data->p[18];
/*data->p[19] =*/ realtype    gbk  ; // data->p[19];
/*data->p[20] =*/ realtype    gk1 ; // data->p[20];
/*data->p[21] =*/ realtype    gkr ; // data->p[21];
/*data->p[22] =*/ realtype    gks ; // data->p[22];
/*data->p[23] =*/ realtype    gcal12 ; // data->p[23];
/*data->p[24] =*/ realtype    gcal13 ; // data->p[24];
/*data->p[25] =*/ realtype    ecal ; // data->p[25];
/*data->p[26] =*/ realtype    kmfca ; // data->p[26];
/*data->p[27] =*/ realtype    alpha_fca ; // data->p[27];
/*data->p[28] =*/ realtype    ecat ; // data->p[28];
/*data->p[29] =*/ realtype    enattxr ; // data->p[29];
/*data->p[30] =*/ realtype    gsus ; // data->p[30];
/*data->p[31] =*/ realtype    inakmax ; // data->p[31];
/*data->p[32] =*/ realtype    kmnap ; // data->p[32];
/*data->p[33] =*/ realtype    kmkp ; // data->p[33];
/*data->p[34] =*/ realtype    kNaCa ; // data->p[34];
/*data->p[35] =*/ realtype    K1ni ; // data->p[35];
/*data->p[36] =*/ realtype    K1no ; // data->p[36];
/*data->p[37] =*/ realtype    K2ni ; // data->p[37];
/*data->p[38] =*/ realtype    K2no ; // data->p[38];
/*data->p[39] =*/ realtype    K3ni ; // data->p[39];
/*data->p[40] =*/ realtype    K3no ; // data->p[40];
/*data->p[41] =*/ realtype    Kci ; // data->p[41];
/*data->p[42] =*/ realtype    Kco ; // data->p[42];
/*data->p[43] =*/ realtype    Kcni ; // data->p[43];
/*data->p[44] =*/ realtype    Qci ; // data->p[44];
/*data->p[45] =*/ realtype    Qco ; // data->p[45];
/*data->p[46] =*/ realtype    Qn ; // data->p[46];
/*data->p[47] =*/ realtype    tdifca ; // data->p[47];
/*data->p[48] =*/ realtype    Prel ; // data->p[48];
/*data->p[49] =*/ realtype    Krel ; // data->p[49];
/*data->p[50] =*/ realtype    nrel ; // data->p[50];
/*data->p[51] =*/ realtype    Kup ; // data->p[51];
/*data->p[52] =*/ realtype    nup ; // data->p[52];
/*data->p[53] =*/ realtype    Ttr ; // data->p[53];
/*data->p[54] =*/ realtype    ConcTC ; // data->p[54];
/*data->p[55] =*/ realtype    ConcTMC ; // data->p[55];
/*data->p[56] =*/ realtype    kfTC ; // data->p[56];
/*data->p[57] =*/ realtype    kfTMC ; // data->p[57];
/*data->p[58] =*/ realtype    kbTC ; // data->p[58];
/*data->p[59] =*/ realtype    kbTMC ; // data->p[59];
/*data->p[60] =*/ realtype    kfTMM ; // data->p[60];
/*data->p[61] =*/ realtype    kbTMM ; // data->p[61];
/*data->p[62] =*/ realtype    ConcCM ; // data->p[62];
/*data->p[63] =*/ realtype    kfCM ; // data->p[63];
/*data->p[64] =*/ realtype    kbCM ; // data->p[64];
/*data->p[65] =*/ realtype    ConcCQ ; // data->p[65];
/*data->p[66] =*/ realtype    kfCQ ; // data->p[66];
/*data->p[67] =*/ realtype    kbCQ ; // data->p[67];
/*data->p[68] =*/ realtype    koca ; // data->p[68];
/*data->p[69] =*/ realtype    kom ; // data->p[69];
/*data->p[70] =*/ realtype    kica ; // data->p[70];
/*data->p[71] =*/ realtype    kim ; // data->p[71];
/*data->p[72] =*/ realtype    eca50sr ; // data->p[72];
/*data->p[73] =*/ realtype    ks ; // data->p[73];
/*data->p[74] =*/ realtype    maxsr ; // data->p[74];
/*data->p[75] =*/ realtype    minsr ; // data->p[75];
/*data->p[76] =*/ realtype    hsrr ; // data->p[76];
/*data->p[77] =*/ realtype    pumpkmf ; // data->p[77];
/*data->p[78] =*/ realtype    pumpkmr ; // data->p[78];
/*data->p[79] =*/ realtype    pumphill ; // data->p[79];
/*data->p[80] =*/ realtype    gna_ttxs ; // data->p[80];
/*data->p[81] =*/ realtype    gna_ttxr ; // data->p[81];
/*data->p[82] =*/ realtype    gcat ; // data->p[82];
/*data->p[83] =*/ realtype    gh ; // data->p[83];
/*data->p[84] =*/ realtype    gto ; // data->p[84];
/*data->p[85] =*/ realtype    Pup ; // data->p[85]; 

/* data->p[86] =  */ realtype  soicr_switch; // , 0.0 , 1 ) // soicr usually is switched off.
/* data->p[87] = */ realtype  soicr_max; // , 25.0 , 1 )
/* data->p[88] =  */ realtype  soicr_thresh; // , 8400.0 , 1 ) // 8.4 mM, i.e. 8400 microM
/* data->p[89] =  */ realtype soicr_slope; // , 1.0 ,  1 ) // 0.001 mM, i.e. 1.0 microM
/* data->p[90] =   */ realtype soicr_tau; // , 125.0 , 1 )

realtype  	   if_vhalf; // 106.8; // 1 ) // ISO
realtype  	   dl13vhalf; // 13.5; // 1 )
realtype  	   fl13vhalf; // 35.0; // 1 )
realtype  	   dl12vhalf; // 3.0; // 1 )
realtype  	   fl12vhalf; // 36.0; // 1 )
realtype  	   ikractvhalf; // 21.173694; // 1 )

realtype cashiftplus;
realtype cashiftminus;
realtype lap_val;

realtype ikr_act_factor;
realtype csqnfree;
realtype open_const;

realtype csqn;


/*data->p[0] =*/         R = data->p[0];
/*data->p[1] =*/         T = data->p[1];
/*data->p[2] =*/         F = data->p[2];
/*data->p[3] =*/         capacitance = data->p[3];
/*data->p[4] =*/         vcell = data->p[4];
/*data->p[5] =*/         l_cell = data->p[5];
/*data->p[6] =*/         r_cell = data->p[6];
/*data->p[7] =*/         vrel = data->p[7];
/*data->p[8] =*/         vsub = data->p[8];
/*data->p[9] =*/         vup  = data->p[9];
/*data->p[10] =*/         vi = data->p[10];
/*data->p[11] =*/         Mgi = data->p[11];
/*data->p[12] =*/         nao = data->p[12];
/*data->p[13] =*/         cao = data->p[13];
/*data->p[14] =*/         ko = data->p[14];
/*data->p[15] =*/         gst = data->p[15];
/*data->p[16] =*/         eist = data->p[16];
/*data->p[17] =*/         gbna = data->p[17];
/*data->p[18] =*/         gbca = data->p[18];
/*data->p[19] =*/         gbk  = data->p[19];
/*data->p[20] =*/         gk1 = data->p[20];
/*data->p[21] =*/         gkr = data->p[21];
/*data->p[22] =*/         gks = data->p[22];
/*data->p[23] =*/         gcal12 = data->p[23];
/*data->p[24] =*/         gcal13 = data->p[24];
/*data->p[25] =*/         ecal = data->p[25];
/*data->p[26] =*/         kmfca = data->p[26];
/*data->p[27] =*/         alpha_fca = data->p[27];
/*data->p[28] =*/         ecat = data->p[28];
/*data->p[29] =*/         enattxr = data->p[29];
/*data->p[30] =*/         gsus = data->p[30];
/*data->p[31] =*/         inakmax = data->p[31];
/*data->p[32] =*/         kmnap = data->p[32];
/*data->p[33] =*/         kmkp = data->p[33];
/*data->p[34] =*/         kNaCa = data->p[34];
/*data->p[35] =*/         K1ni = data->p[35];
/*data->p[36] =*/         K1no = data->p[36];
/*data->p[37] =*/         K2ni = data->p[37];
/*data->p[38] =*/         K2no = data->p[38];
/*data->p[39] =*/         K3ni = data->p[39];
/*data->p[40] =*/         K3no = data->p[40];
/*data->p[41] =*/         Kci = data->p[41];
/*data->p[42] =*/         Kco = data->p[42];
/*data->p[43] =*/         Kcni = data->p[43];
/*data->p[44] =*/         Qci = data->p[44];
/*data->p[45] =*/         Qco = data->p[45];
/*data->p[46] =*/         Qn = data->p[46];
/*data->p[47] =*/         tdifca = data->p[47];
/*data->p[48] =*/         Prel = data->p[48];
/*data->p[49] =*/         Krel = data->p[49];
/*data->p[50] =*/         nrel = data->p[50];
/*data->p[51] =*/         Kup = data->p[51];
/*data->p[52] =*/         nup = data->p[52];
/*data->p[53] =*/         Ttr = data->p[53];
/*data->p[54] =*/         ConcTC = data->p[54];
/*data->p[55] =*/         ConcTMC = data->p[55];
/*data->p[56] =*/         kfTC = data->p[56];
/*data->p[57] =*/         kfTMC = data->p[57];
/*data->p[58] =*/         kbTC = data->p[58];
/*data->p[59] =*/         kbTMC = data->p[59];
/*data->p[60] =*/         kfTMM = data->p[60];
/*data->p[61] =*/         kbTMM = data->p[61];
/*data->p[62] =*/         ConcCM = data->p[62];
/*data->p[63] =*/         kfCM = data->p[63];
/*data->p[64] =*/         kbCM = data->p[64];
/*data->p[65] =*/         ConcCQ = data->p[65];
/*data->p[66] =*/         kfCQ = data->p[66];
/*data->p[67] =*/         kbCQ = data->p[67];
/*data->p[68] =*/         koca = data->p[68];
/*data->p[69] =*/         kom = data->p[69];
/*data->p[70] =*/         kica = data->p[70];
/*data->p[71] =*/         kim = data->p[71];
/*data->p[72] =*/         eca50sr = data->p[72];
/*data->p[73] =*/         ks = data->p[73];
/*data->p[74] =*/         maxsr = data->p[74];
/*data->p[75] =*/         minsr = data->p[75];
/*data->p[76] =*/         hsrr = data->p[76];
/*data->p[77] =*/         pumpkmf = data->p[77];
/*data->p[78] =*/         pumpkmr = data->p[78];
/*data->p[79] =*/         pumphill = data->p[79];
/*data->p[80] =*/         gna_ttxs = data->p[80];
/*data->p[81] =*/         gna_ttxr = data->p[81];
/*data->p[82] =*/         gcat = data->p[82];
/*data->p[83] =*/         gh = data->p[83];
/*data->p[84] =*/         gto = data->p[84];
/*data->p[85] =*/         Pup = data->p[85]; 

soicr_switch = data->p[86]; // , 0.0 , 1 ) // soicr usually is switched off.
soicr_max  = data->p[87]; // , 25.0 , 1 )
soicr_thresh  = data->p[88]; // , 8400.0 , 1 ) // 8.4 mM, i.e. 8400 microM
soicr_slope  = data->p[89]; // , 1.0 ,  1 ) // 0.001 mM, i.e. 1.0 microM
soicr_tau  = data->p[90]; // , 125.0 , 1 )

/* realtype */  	   if_vhalf = data->p[91]; // 106.8 = data->p[]; // 1 ) // ISO
/* realtype */  	   dl13vhalf = data->p[92]; // 13.5 = data->p[]; // 1 )
/* realtype */  	   fl13vhalf = data->p[93]; // 35.0 = data->p[]; // 1 )
/* realtype */  	   dl12vhalf = data->p[94]; // 3.0 = data->p[]; // 1 )
/* realtype */  	   fl12vhalf = data->p[95]; // 36.0 = data->p[]; // 1 )
/* realtype */  	   ikractvhalf = data->p[96]; // 21.173694 = data->p[]; // 1 )
			   cashiftplus = data->p[97];
			   lap_val = data->lap_val;
			   cashiftminus = data->p[98];

			   ikr_act_factor = data->p[99];
	
			   csqnfree = data->p[100];
			   open_const = data->p[101];

/*
the formulas and the RHS for KLYZMOUSE
*/

/* Model variables */

// double v;
double sstime;
int  current_trace_counter = 0, output_counter = 0;

double vclamp;
double total_current;

/*
Ion concentrations and reversal potentials
*/

//double cai;             // mM
//double casub;           // mM 

double ena,eca,ek,eks;

/*
Ist variables.
*/

double ist;
double qa, qi, tauqa,tauqi,alphaqa,betaqa,alphaqi,betaqi;
//double dst; // ist gating variables
//double fst;

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
//double dt;
//double ft;
double dt_inf, tau_dt, ft_inf, tau_ft;

/*
IKr variables
*/

double ikr_act_inf,tau_ikr_act;
double ikr;
double ikr_inact_inf,tau_ikr_inact,tau_ikr_inact2;
// ,ikr_inact,ikr_inact2;

/*
IKs variables
*/

// double iks_act
double iks_act_inf,tau_iks_act;
double iks;

/*
ICaL 1.2 and 1.3 parameters.
*/

double alpha_dl, beta_dl, tau_dl, alpha_fl, beta_fl, tau_fl;
// double fl12, dl12, 
double dl12_inf, fl12_inf,ical12;
// double fl13, dl13, 
double dl13_inf, fl13_inf,ical13;

// double fca; 
double fca_inf ;
double taufca ;

/*
INa Nav1.1 (TTXS) and Nav1.5 (TTXR) variables,
*/

double ina_ttxr, ina_ttxs; // ina: Nav 1.1 and Nav 1.5
double m3_inf_ttxr, h_inf_ttxr;
double m3_inf_ttxs, h_inf_ttxs;
double m_inf_ttxr,m_inf_ttxs, hs,hsr;
//double m_ttxr,h_ttxr,j_ttxr;
//double m_ttxs,h_ttxs,j_ttxs;
double tau_m,tau_h,tau_j,tau_mr,tau_hr,tau_jr;
double delta_ttxr_m,delta_ttxr_h,delta_ttxr_j;
double fna;

/*
If variables
*/

double ih, ihk, ihna, ih_1, ih_2, ih_4, /* y_1_2, */ y_4, tau_y_1_2, tau_y_4, y_inf, y_1_2_inf, y_4_inf,ih_1_k, ih_2_k, ih_4_k, ih_1_na, ih_2_na, ih_4_na;

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

double r_inf, /* r,*/ tau_r, isus;

/*
Ito variables
*/

// double q,
double q_inf,tau_q;
double ito;

/*
Ionic homeostasis
*/

double ca_flux;

/*
Calcium diffusion
*/

double Jcadif;

double Jrel, Jup, Jtr;
//,carel,caup;

//double dFtc; // ,Ftc;
//double dFtmc; // ,Ftmc,Ftmm;
//double dFtmm; // ,Ftmm;
//double dFcms; // ,Fcms;
//double dFcmi; // ,Fcmi;
//double dFcq; // ,Fcq;

//double nai=8.0;         // mM 
//double ki=140;
double FRT, RTF;

/*
RyR markov chain
*/

// double resting, open, inactivated, resting_inactivated;
double kcasr,kosrca,kisrca;

double nai_tot, ki_tot;
double dvdt;

double soicr_ss;

#define v  	(Y[0])
#define dst  	(Y[1])
#define fst  	(Y[2])
#define dt  	(Y[3])
#define ft  	(Y[4])
#define ikr_act  (Y[5])
#define ikr_inact  (Y[6])
#define ikr_inact2  (Y[7])
#define iks_act  (Y[8])
#define fl12  	(Y[9])
#define dl12  	(Y[10])
#define fl13  	(Y[11])
#define dl13  	(Y[12])
#define r  	(Y[13])
#define m_ttxr  (Y[14])
#define h_ttxr  (Y[15])
#define j_ttxr  (Y[16])
#define m_ttxs  (Y[17])
#define h_ttxs  (Y[18])
#define j_ttxs  (Y[19])
#define y_1_2  	(Y[20])
#define carel  	(Y[21])
#define caup  	(Y[22])
#define casub  	(Y[23])
#define Ftc  	(Y[24])
#define Ftmc  	(Y[25])
#define Ftmm  	(Y[26])
#define Fcms  	(Y[27])
#define Fcmi  	(Y[28])
#define Fcq  	(Y[29])
#define cai  	(Y[30])
#define q  	(Y[31])
#define fca  	(Y[32])
#define nai  	(Y[33])
#define ki  	(Y[34])
#define resting  (Y[35])
#define open  	(Y[36])
#define inactivated (Y[37])
#define resting_inactivated (Y[38])
#define soicr (Y[39])

double ddst, dfst, dddt, dft, dikr_act, dikr_inact, dikr_inact2, diks_act, dfl12, ddl12, dfl13, ddl13, dr, dm_ttxr,
	dh_ttxr, dj_ttxr, dm_ttxs, dh_ttxs, dj_ttxs, dy_1_2, dcarel, dcaup, dcasub, dFtc, dFtmc, dFtmm, dFcms, 
	dFcmi, dFcq, dcai, dq, dfca, dnai, dki, dresting, dopen, dinactivated, dresting_inactivated;
double dsoicr;

FRT = F/(R*T);

ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);

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

ddst = (qa-dst)/tauqa;
dfst = (qi-fst)/tauqi;

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
dddt = (dt_inf - dt)/tau_dt;

tau_ft = 1.0/(0.0153*exp(-(v+61.7)/83.3)+0.015*exp((v+61.7)/15.38)); // Kurata/Maltsev.
ft_inf = 1.0/(1.0+exp((v + 61.7)/5.6)); // rabbit SS, since the Mangoni data is 10 mV low, and very steep.
dft = (ft_inf - ft)/tau_ft;

icat = gcat*ft*dt*(v - ecat);           // nA

/*Ikr*************************************************************************/

ikr_act_inf = 1.0/(1.0 + exp(-(v+ikractvhalf)/9.757086)); // sacred // ISO shifts it to the left
tau_ikr_act = ikr_act_factor*0.699821/(0.003596*exp((v)/15.339290) + 0.000177*exp(-(v)/25.868423)); // this fits the Q10 variation.
dikr_act = (ikr_act_inf-ikr_act)/tau_ikr_act;
     
/* SS inactivation */

ikr_inact_inf = 1.0/(1.0 + exp((v+20.758474-4.0)/(19.0)));
tau_ikr_inact = 0.2+0.9*1.0/(0.1*exp(v/54.645)+0.656*exp(v/106.157));
dikr_inact = (ikr_inact_inf - ikr_inact)/tau_ikr_inact;
dikr_inact2 = 0.0; // this variable is not used - developmental (2010) artefact.
ikr = gkr*ikr_act*ikr_inact*(v - ek);

/**IKs********************************************************************/

iks_act_inf = 1.0/(1.0 + exp(-(v-20.876040)/11.852723));
tau_iks_act =  1000.0/(13.097938/(1.0 + exp(-(v-48.910584)/10.630272)) + exp(-(v)/35.316539)); // Zhang model. The other models are similar.
diks_act = (iks_act_inf - iks_act)/tau_iks_act; // Ding/Maatsura
iks = gks*iks_act*iks_act*(v - eks);

/*ICaL*******************************************************************/

/*
realtype  	   dl13vhalf; // 13.5; // 1 )
realtype  	   fl13vhalf; // 35.0; // 1 )
realtype  	   dl12vhalf; // 3.0; // 1 )
realtype  	   fl12vhalf; // 36.0; // 1 )
*/

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

dl13_inf = 1.0/(1+exp(-(v+dl13vhalf)/6.0));
fl13_inf = 1.0/(1+exp((v+fl13vhalf)/7.3));
tau_fl = (7.4 + 45.77*exp(-0.5*(v+28.1)*(v+28.1)/(11*11)));

ddl13 = (dl13_inf - dl13)/tau_dl;
dfl13 = (fl13_inf - fl13)/tau_fl;

/* Cav 1.2 */

dl12_inf = 1.0/(1+exp(-(v+dl12vhalf)/5.0)); // Mangoni 2006: according to his email
fl12_inf = 1.0/(1+exp((v+fl12vhalf)/4.6)); // Mangoni 2006 35.9

ddl12 = (dl12_inf - dl12)/tau_dl;
dfl12 = (fl12_inf - fl12)/tau_fl;

/* fca */

fca_inf = kmfca/(kmfca+casub);
taufca = fca_inf/alpha_fca;
dfca = (fca_inf - fca)/taufca;

 ical12 = gcal12*fl12*dl12*fca*(v-ecal);

 ical13 = gcal13*fl13*dl13*fca*(v-ecal);
    
/**INa**********************************************************************/

fna = (9.52e-02*exp(-6.3e-2*(v+34.4))/(1+1.66*exp(-0.225*(v+63.7))))+8.69e-2; // 1-fna = proportion of fast inactivation, fna = contribution of slow

m3_inf_ttxr = 1.0/(1.0 + exp(-(v+45.213705)/7.219547));
h_inf_ttxr = 1.0/(1.0 + exp(-(v+62.578120 )/(-6.084036)));
m3_inf_ttxs = 1.0/(1.0 + exp(-(v+36.097331-5.0)/5.0)); // Mangoni takes -29.
h_inf_ttxs = 1.0/(1.0 + exp((v+56.0)/3.0));

m_inf_ttxr = pow(m3_inf_ttxr,0.333);
m_inf_ttxs = pow(m3_inf_ttxs,0.333);

tau_m = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_h = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_j = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);

dm_ttxs = (m_inf_ttxs - m_ttxs)/tau_m;
dh_ttxs = (h_inf_ttxs - h_ttxs)/tau_h;
dj_ttxs = (h_inf_ttxs - j_ttxs)/tau_j;
 
hs = (1.0-fna)*h_ttxs+fna*j_ttxs;

tau_mr = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_hr = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_jr = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);

dm_ttxr = (m_inf_ttxr - m_ttxr)/tau_mr;
dh_ttxr = (h_inf_ttxr - h_ttxr)/tau_hr;
dj_ttxr = (h_inf_ttxr - j_ttxr)/tau_jr;

hsr = (1.0-fna)*h_ttxr+fna*j_ttxr;

if(fabs(v)>0.005)
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*(F*F/(R*T))*((exp((v-ena)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*F*((exp((v-ena)*F/(R*T))-1.0));


if(fabs(v)>0.005)
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*(F*F/(R*T))*((exp((v-enattxr)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*F*((exp((v-enattxr)*F/(R*T))-1.0));


/**If**************************************************************************/

y_inf = 1.0/(1.0 + exp((v+if_vhalf)/16.3));
tau_y_1_2 = 1.5049/(exp(-(v+590.3)*0.01094)+ exp((v-85.1)/17.2));
dy_1_2 = (y_inf - y_1_2)/tau_y_1_2;

ihk  = 0.6167*gh*y_1_2*(v - ek);
ihna = 0.3833*gh*y_1_2*(v - ena);

ih = (ihk + ihna);

/*Ito*************************************************************************/

q_inf = 1.0/(1.0+exp((v+49.0)/13.0));
tau_q = (6.06 + 39.102/(0.57*exp(-0.08*(v+44.0))+0.065*exp(0.1*(v+45.93))))/0.67; 
dq = ((q_inf-q)/tau_q);
r_inf = 1.0/(1.0+exp(-(v-19.3)/15.0));
tau_r = (2.75+14.40516/(1.037*exp(0.09*(v+30.61))+0.369*exp(-0.12*(v+23.84))))/0.303;
dr = (r_inf-r)/tau_r;
ito = gto*q*r*(v-ek);
 
/*Isus***********************************************************************/

isus = gsus*r*(v-ek);

/*Inak***********************************************************************/

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

 kosrca = cashiftplus*koca/kcasr;

 kisrca = cashiftminus*kica*kcasr;

dresting = kim*resting_inactivated - kisrca*casub*resting - kosrca*casub*casub*resting + kom*open;
dopen    = kosrca*casub*casub*resting - kom*open - kisrca*casub*open + kim*inactivated;
dinactivated = kisrca*casub*open - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated;
dresting_inactivated = kom*inactivated - kosrca*casub*casub*resting_inactivated - kim*resting_inactivated + kisrca*casub*resting;

dsoicr = 0.0; // (soicr_ss - soicr)/soicr_tau; // do this in the mousesk.c, not f()

Jrel = (ks+soicr_switch*soicr)*open*(carel - casub); // original formula.

// Jrel = ks*(open+soicr_switch*soicr)*(carel - casub); // this is for the soicr story so far, modified.

// if(csqnfree>0)
// Jrel = ks*open_const*(carel - casub); // this is more effective, at least in the control case.

/*
basal: this number is between 7.5 and 10.
ISO without messabout: this number is between 8.5 and 10.
In the next few lines, csqnfree is used to denote the csqnthreshold,
as well as the switching on of the csqn calculation.
*/
/*
if(csqnfree>0.0){

	csqn = ConcCQ*(1.0 - Fcq);
	if(csqnfree<csqn)
	Jrel = ks*open_const*(carel - casub);
	else
	Jrel = ks*(open+soicr_switch*soicr)*(carel - casub);

}
*/
/*
Carel based release.
This is the most logical among what I have read so far.
The literature says that under ISO, the Carel is increased.
With ISO and the ca2+ release disease, the krel is increased.
start by modelling it simplistically as an if statement.
I am using csqnfree to represent carel_threshold
*/
/*
	if(carel>csqnfree)
	Jrel = ks*open_const*(carel - casub);
	else
	Jrel = ks*(carel - casub);
*/

Jup = Pup*(pow(cai/pumpkmf,pumphill) - pow(caup/pumpkmr,pumphill))/(1.0 + pow(cai/pumpkmf,pumphill) + pow(caup/pumpkmr,pumphill));

Jtr  = (caup - carel)/Ttr;

// Ca buffering flux // these are the derivatives of the F's
dFtc  = kfTC*cai*(1.0-Ftc)-kbTC*Ftc;
dFtmc = kfTMC*cai*(1.0-Ftmc-Ftmm)-kbTMC*Ftmc;
dFtmm = kfTMM*Mgi*(1.0-Ftmc-Ftmm)-kbTMM*Ftmm;
dFcms = kfCM*casub*(1.0-Fcms)-kbCM*Fcms;
dFcmi = kfCM*cai*(1.0-Fcmi)-kbCM*Fcmi;
dFcq  = kfCQ*carel*(1.0-Fcq)-kbCQ*Fcq;

dcasub = ((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);
dcai = ((Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmi + ConcTC*dFtc + ConcTMC*dFtmc)); 

dcarel = (Jtr - Jrel - ConcCQ*dFcq);
dcaup = (Jup-Jtr*vrel/vup);

/***************************************************************************/

total_current = ih+ina_ttxr+ina_ttxs+ical12+ical13+iks+ikr+ik1+ist+ib+icat+inak+isus+inaca+ito;
dvdt = /* lap_val/capacitance + */ (- total_current)/capacitance;
 
ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);
 
nai_tot = ihna+ina_ttxr+ina_ttxs+3.0*inak+3.0*inaca+ist+ibna;
ki_tot = ihk+iks+ikr+ik1+ibk-2.0*inak+isus+ito;

dnai = -(nai_tot)/(F*vi);
dki = -(ki_tot)/(F*vi);

dY[0] = dvdt;
dY[1] = ddst;
dY[2] = dfst;
dY[3] = dddt;
dY[4] = dft;
dY[5] = dikr_act;
dY[6] = dikr_inact;
dY[7] = dikr_inact2;
dY[8] = diks_act;
dY[9] = dfl12;
dY[10] = ddl12;
dY[11] = dfl13;
dY[12] = ddl13;
dY[13] = dr;
dY[14] = dm_ttxr;
dY[15] = dh_ttxr;
dY[16] = dj_ttxr;
dY[17] = dm_ttxs;
dY[18] = dh_ttxs;
dY[19] = dj_ttxs;
dY[20] = dy_1_2;
dY[21] = dcarel;
dY[22] = dcaup;
dY[23] = dcasub;
dY[24] = dFtc;
dY[25] = dFtmc;
dY[26] = dFtmm;
dY[27] = dFcms;
dY[28] = dFcmi;
dY[29] = dFcq;
dY[30] = dcai;
dY[31] = dq;
dY[32] = dfca;
dY[33] = dnai;
dY[34] = dki;
dY[35] = dresting;
dY[36] = dopen;
dY[37] = dinactivated;
dY[38] = dresting_inactivated;

/***********************************************************************/

 for(i=0;i<NEQS;i++)
  Ith(ydot,i+1) = dY[i];

 return(0);
}

