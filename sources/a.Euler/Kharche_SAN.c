/*
July 2011.
A mathematical model of action potentials of mouse sinoatrial node cells with molecular bases
Sanjay Kharche, Jian Yu, Ming Lei, and Henggui Zhang.
AJP - Heart September 2011 vol. 301 no. 3 H945-H963

LINK TO PAPER:

http://ajpheart.physiology.org/content/301/3/H945.long
*/

#include <stdio.h>
#include <math.h>
#include<stdlib.h>

#define R 8.314472
#define T 310.5   
#define F 96.4845

#define number_of_apds 100

/***************************************************************************************/

#define ddt 0.00025
#define capacitance 0.025
#define vcell 3.0
#define l_cell 66.3767257
#define r_cell 3.792956
#define vrel 0.0036
#define vsub 0.03328117
#define vup  0.0348
#define vi 1.34671883
#define Mgi 2.5
#define nao 140.0
#define cao 1.8
#define ko 5.4
#define gst 0.00006
#define eist 17.0
#define gbna 0.0001215
#define gbca 0.000015
#define gbk  0.0000025
#define	gk1 0.229*0.0039228*0.9
#define gks 0.000299
#define ecal 47.0
#define kmfca 0.00035
#define alpha_fca 0.021
#define all_ica_multiplier 1.0
#define ecat 45.0
#define enattxr  41.5761
#define multiplier2 1.0
#define gsus 0.00039060
#define inakmax_multiplier 1.85
#define inakmax inakmax_multiplier*0.077
#define kmnap 14.0
#define kmkp 1.4
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
#define tdifca 0.04
#define Prel 2.5
#define Krel 0.0015
#define nrel 2.0
#define Kup 0.0006
#define nup 1.0
#define Ttr 40.0
#define ConcTC 0.031
#define ConcTMC 0.062
#define kfTC 88.8
#define kfTMC 237.7
#define kbTC 0.446
#define kbTMC 0.00751
#define kfTMM 2.277
#define kbTMM 0.751
#define ConcCM 0.045
#define kfCM 237.7
#define kbCM 0.542
#define ConcCQ 10.0
#define kfCQ 0.534
#define kbCQ 0.445
#define koca 10.0
#define kom 0.06
#define kica 0.5
#define kim 0.005
#define eca50sr 0.45
#define maxsr 15.0
#define minsr 1.0
#define hsrr 2.5
#define pumphill 2.0

double v, vnew,sstime;
int  current_trace_counter = 0, output_counter = 0;
double vclamp;
double total_current;
double cai; 
double casub;
double ena,eca,ek,eks;
double ist;
double qa, qi, tauqa,tauqi,alphaqa,betaqa,alphaqi,betaqi;
double dst;
double fst;
double ibca,ibna,ibk,ib;
double ik1 ;
double xk1inf;
double icat;
double dt;
double ft;
double dt_inf, tau_dt, ft_inf, tau_ft;
double ikr_act, ikr_act_inf,tau_ikr_act;
double ikr;
double ikr_inact_inf,tau_ikr_inact,tau_ikr_inact2,ikr_inact,ikr_inact2;
double iks_act, iks_act_inf,tau_iks_act;
double iks;
double alpha_dl, beta_dl, tau_dl, alpha_fl, beta_fl, tau_fl;
double fl12, dl12, dl12_inf, fl12_inf,ical12;
double fl13, dl13, dl13_inf, fl13_inf,ical13;
double fca; 
double fca_inf;
double taufca;
double ina_ttxr, ina_ttxs;
double m3_inf_ttxr, h_inf_ttxr;
double m3_inf_ttxs, h_inf_ttxs;
double m_inf_ttxr,m_inf_ttxs, hs,hsr;
double m_ttxr,h_ttxr,j_ttxr;
double m_ttxs,h_ttxs,j_ttxs;
double tau_m,tau_h,tau_j,tau_mr,tau_hr,tau_jr;
double delta_ttxr_m,delta_ttxr_h,delta_ttxr_j;
double fna;
double ih, ihk, ihna, ih_1, ih_2, ih_4, y_1_2, y_4, tau_y_1_2, tau_y_4, y_inf, y_1_2_inf, y_4_inf,ih_1_k, ih_2_k, ih_4_k, ih_1_na, ih_2_na, ih_4_na;
double inak;
double inaca;
double di,doo,k43,k12,k14,k41,k34,k21,k23,k32,x1,x2,x3,x4;
double r_inf, r, tau_r, isus;
double q,q_inf,tau_q;
double ito;
double ca_flux;
double Jcadif;
double Jrel, Jup, Jtr,carel,caup;
double variation_multiply = 1.0;
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
double nai=8.0;  
double ki=140;
double FRT, RTF;
double resting, open, inactivated, resting_inactivated;
double kcasr,kosrca,kisrca;
double pumpvmax = 0.04;
double pumpkmf = 0.000246;
double pumpkmr = 3.29;
double pumpfactor = 1.0;
double Pup;
double Pup2;
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

int main()
{  
FILE *outputcurrents,*parameters, *random;

char *str;
int start_output = 0;
double total_time = 0.0;
int temp_counter = 0;

FRT = F/(R*T);
param_counter = 0;
current_trace_counter = 0;
outputcurrents = fopen("mousesan.dat","w+");

/*
Initial conditions
*/

v = -64.5216286940;
dst = 0.6246780312;
fst = 0.4537033169;
dt = 0.0016256324;
ft = 0.4264459666;
ikr_act = 0.4043600437;
ikr_inact = 0.9250035423;
ikr_inact2 = 0.1875749806;
iks_act = 0.0127086259;
fl12 = 0.9968141226;
dl12 = 0.0000045583;
fl13 = 0.9809298233;
dl13 = 0.0002036671;
r = 0.0046263658;
m_ttxr = 0.4014088304;
h_ttxr = 0.2724817537;
j_ttxr = 0.0249208708;
m_ttxs = 0.1079085266;
h_ttxs = 0.4500098710;
j_ttxs = 0.0268486392;
y_1_2 = 0.0279984462;
y_4 = 0.0137659036;
carel = 0.1187281829;
caup = 1.5768287365;
casub = 0.0000560497;
Ftc = 0.0063427103;
Ftmc = 0.1296677919;
Ftmm = 0.7688656371;
Fcms = 0.0242054739;
Fcmi = 0.0138533048;
Fcq = 0.1203184861;
cai = 0.0000319121;
q = 0.6107148187;
fca = 0.7649576191;
nai = 8.1179761505;
ki = 139.8854603066;
resting = 0.7720290515;
open = 0.0000000760;
inactivated = 0.0000000213;
resting_inactivated = 0.2162168926;

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
Pup2 = 0.04;
pumpkmf = 0.00008;
pumpkmr = 4.5;

for(param_counter=0;param_counter<number_of_apds;param_counter++){
     min_potential[param_counter] = 100000.0;
     max_potential[param_counter] = -10000.0;
     dvdtmax[param_counter]       = -10000.0;
     ddr[param_counter]           = -10000.0;
     top[param_counter]           =  10000.0;
     top_slope[param_counter]     = -10000.0;
     apd50[param_counter]         = -10000.0;
     apd90[param_counter]         = -10000.0;
     cycle_length[param_counter]  = -10000.0;
    }

ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);

param_counter = 0;
dvdt = 1000.0;
dvdtold = 500.0;
start_output = 0;
total_time = 5000.0;

for(sstime=0.0;sstime<=total_time;sstime = sstime + ddt){

/*Ist********************************************************************/

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

tau_dt = 1.0/(1.068*exp((v + 26.3)/30.0) + 1.068*exp(-(v + 26.3)/30.0));
dt_inf = 1.0/(1.0+exp(-(v + 26.0)/6.0));
dt = dt + ddt*((dt_inf - dt)/tau_dt);
tau_ft = 1.0/(0.0153*exp(-(v+61.7)/83.3)+0.015*exp((v+61.7)/15.38));
ft_inf = 1.0/(1.0+exp((v + 61.7)/5.6));
ft = ft + ddt*((ft_inf - ft)/tau_ft);
icat = gcat*ft*dt*(v - ecat);          

/*Ikr********************************************************************/

ikr_act_inf = 1.0/(1.0 + exp(-(v+21.173694)/9.757086));
tau_ikr_act = 0.699821/(0.003596*exp((v)/15.339290) + 0.000177*exp(-(v)/25.868423));
ikr_act = ikr_act + ddt*(ikr_act_inf-ikr_act)/tau_ikr_act;     
ikr_inact_inf = 1.0/(1.0 + exp((v+20.758474-4.0)/(19.0)));
tau_ikr_inact = 0.2+0.9*1.0/(0.1*exp(v/54.645)+0.656*exp(v/106.157));
ikr_inact = ikr_inact + ddt*(ikr_inact_inf - ikr_inact)/tau_ikr_inact;
ikr = gkr*ikr_act*ikr_inact*(v - ek);

/**IKs********************************************************************/

iks_act_inf = 1.0/(1.0 + exp(-(v-20.876040)/11.852723));
tau_iks_act =  1000.0/(13.097938/(1.0 + exp(-(v-48.910584)/10.630272)) + exp(-(v)/35.316539));
iks_act = iks_act + ddt*(iks_act_inf - iks_act)/tau_iks_act;
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
dl13_inf = 1.0/(1+exp(-(v+13.5)/6.0));
fl13_inf = 1.0/(1+exp((v+35.0)/7.3));
tau_fl = (7.4 + 45.77*exp(-0.5*(v+28.1)*(v+28.1)/(11*11)));
dl13 = dl13 + ddt*(dl13_inf - dl13)/tau_dl;
fl13 = fl13 + ddt*(fl13_inf - fl13)/tau_fl;
dl12_inf = 1.0/(1+exp(-(v+3.0)/5.0));
fl12_inf = 1.0/(1+exp((v+36.0)/4.6));
dl12 = dl12 + ddt*(dl12_inf - dl12)/tau_dl;
fl12 = fl12 + ddt*(fl12_inf - fl12)/tau_fl;
fca_inf = kmfca/(kmfca+casub);
taufca = fca_inf/alpha_fca;
fca = fca + ddt*(fca_inf - fca)/taufca;
ical12 = gcal12*fl12*dl12*fca*(v-ecal);
ical13 = gcal13*fl13*dl13*fca*(v-ecal);
   
/**INa**********************************************************************/

fna = (9.52e-02*exp(-6.3e-2*(v+34.4))/(1+1.66*exp(-0.225*(v+63.7))))+8.69e-2; 
m3_inf_ttxr = 1.0/(1.0 + exp(-(v+45.213705)/7.219547));
h_inf_ttxr = 1.0/(1.0 + exp(-(v+62.578120 )/(-6.084036)));
m3_inf_ttxs = 1.0/(1.0 + exp(-(v+36.097331-5.0)/5.0));
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
if(fabs(v)>0.005)
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*(F*F/(R*T))*((exp((v-ena)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*F*((exp((v-ena)*F/(R*T))-1.0));
if(fabs(v)>0.005)
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*(F*F/(R*T))*((exp((v-enattxr)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*F*((exp((v-enattxr)*F/(R*T))-1.0));

/**If**************************************************************************/

y_inf = 1.0/(1.0 + exp((v+vhalf_gh)/16.3));
tau_y_1_2 = 1.5049/(exp(-(v+590.3)*0.01094)+ exp((v-85.1)/17.2));
y_1_2 = y_1_2 + ddt*(y_inf - y_1_2)/tau_y_1_2;
ihk  = 0.6167*gh*y_1_2*(v - ek);
ihna = 0.3833*gh*y_1_2*(v - ena);
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
inak = inakmax*(pow(ko,1.2)/(pow(kmkp,1.2)+pow(ko,1.2)))*(pow(nai,1.3)/(pow(kmnap,1.3)+pow(nai,1.3)))/(1.0+exp(-(v-ena+120.0)/30.0));

/****iNaCa*******************************************************************/

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
ca_flux = (ical12+ical13+icat-2.0*inaca+ibca)/(2.0*F);
Jcadif = (casub - cai)/tdifca;
kcasr = maxsr - (maxsr - minsr)/(1.0 + pow(eca50sr/carel,hsrr));
kosrca = koca/kcasr;
kisrca = kica*kcasr;
resting = resting + ddt*(kim*resting_inactivated - kisrca*casub*resting - kosrca*casub*casub*resting + kom*open);
open    = open    + ddt*(kosrca*casub*casub*resting - kom*open - kisrca*casub*open + kim*inactivated);
inactivated = inactivated + ddt*(kisrca*casub*open - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated);
resting_inactivated = resting_inactivated + ddt*(kom*inactivated - kosrca*casub*casub*resting_inactivated - kim*resting_inactivated + kisrca*casub*resting);
Jrel = ks*open*(carel - casub);
Jup = Pup*(pow(cai/pumpkmf,pumphill) - pow(caup/pumpkmr,pumphill))/(1.0 + pow(cai/pumpkmf,pumphill) + pow(caup/pumpkmr,pumphill));
Jtr  = (caup - carel)/Ttr;
dFtc  = kfTC*cai*(1.0-Ftc)-kbTC*Ftc;
dFtmc = kfTMC*cai*(1.0-Ftmc-Ftmm)-kbTMC*Ftmc;
dFtmm = kfTMM*Mgi*(1.0-Ftmc-Ftmm)-kbTMM*Ftmm;
dFcms = kfCM*casub*(1.0-Fcms)-kbCM*Fcms;
dFcmi = kfCM*cai*(1.0-Fcmi)-kbCM*Fcmi;
dFcq  = kfCQ*carel*(1.0-Fcq)-kbCQ*Fcq;
Ftc = Ftc + ddt*dFtc;        
Ftmc = Ftmc + ddt*dFtmc;     
Ftmm = Ftmm + ddt*dFtmm;     
Fcms = Fcms + ddt*dFcms;     
Fcmi = Fcmi + ddt*dFcmi;     
Fcq = Fcq + ddt*dFcq;        
casub = casub + ddt*((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);
cai = cai + ddt*((Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmi + ConcTC*dFtc + ConcTMC*dFtmc)); 
carel = carel + ddt*(Jtr - Jrel - ConcCQ*dFcq);
caup = caup + ddt*(Jup-Jtr*vrel/vup);
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
nai = nai - ddt*(nai_tot)/(F*vi);
ki = ki - ddt*(ki_tot)/(F*vi);

if(dvdt>=0.0&&dvdtold<0.0){
 min_potential[param_counter] = v;
 nai_min = nai;
 ki_min = ki;
 tmin_potential[param_counter] = sstime;
 start_output = 1;
}
if(dvdt>dvdtmax[param_counter]&&start_output>0){
 dvdtmax[param_counter] = dvdt;
 apd_start[param_counter] = sstime;
 vdvdtmax[param_counter] = v;	
}
if(dvdtold>0.0&&dvdt<=0.0){
 max_potential[param_counter] = v;
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
  printf("%10.5f\t",min_potential[param_counter]);
  printf("%10.5f\t",max_potential[param_counter]);
  printf("%10.5f\t",dvdtmax[param_counter]);
  printf("%10.5f\t",apd50[param_counter]);
  printf("%10.5f\t",apd90[param_counter]);
  printf("%10.5f\t",cycle_length[param_counter]);
  printf("%10.5f\t",ddr[param_counter]);
  printf("%10.5f\t",top[param_counter]);
  printf("\n");
  param_counter++;
 }
}
v = vnew;

if(current_trace_counter%10==0){
fprintf(outputcurrents,"%10.10f\t%10.10f\t",sstime,v);
fprintf(outputcurrents,"\n");
}

current_trace_counter++;

} // end of time loop

fclose(outputcurrents);

return 0;
}
