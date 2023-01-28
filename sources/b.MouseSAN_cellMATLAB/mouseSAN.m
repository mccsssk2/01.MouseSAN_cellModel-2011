function dumi = mouseSAN()
%
% Sanjay R. Kharche.
% 30 July 2020.
% Tested by Kapiraj Chandrabalan (kapiraj.chandrabalan@lhsc.on.ca).
% This is the mouse sinoatrial node cell model (Am J Physiol Heart Circ Physiol. 2011 Sep; 301(3): H945–H963).
%
% Inputs: none.
%
DELTAT = 0.00025; % units of time: ms.
totalTime = 300.0; % ms.
% initial conditions of state variables.
y(1) 	= -64.5216286940; % v = -64.5216286940
y(2) 	= 0.6246780312; % dst = 0.6246780312
y(3) 	= 0.4537033169; % fst = 0.4537033169
y(4) 	= 0.0016256324; % dt = 0.0016256324
y(5)	= 0.4264459666; % ft = 0.4264459666
y(6) 	= 0.4043600437; % ikr_act = 0.4043600437;
y(7) 	= 0.9250035423; % ikr_inact = 0.9250035423;
y(8) 	= 0.1875749806; % ikr_inact2 = 0.1875749806;
y(9) 	= 0.0127086259; % iks_act = 0.0127086259;
y(10) 	= 0.9968141226; % fl12 = 0.9968141226;
y(11) 	= 0.0000045583; % dl12 = 0.0000045583;
y(12) 	= 0.9809298233; % fl13 = 0.9809298233;
y(13) 	= 0.0002036671; % dl13 = 0.0002036671;
y(14) 	= 0.0046263658; % r = 0.0046263658;
y(15) 	= 0.4014088304; % m_ttxr = 0.4014088304
y(16)	= 0.2724817537; % h_ttxr = 0.2724817537
y(17) 	= 0.0249208708; % j_ttxr = 0.0249208708
y(18) 	= 0.1079085266; % m_ttxs = 0.1079085266
y(19) 	= 0.4500098710; % h_ttxs = 0.4500098710
y(20) 	= 0.0268486392; % j_ttxs = 0.0268486392
y(21) 	= 0.0279984462; % y_1_2 = 0.0279984462
y(22) 	= 0.0137659036; % y_4 = 0.0137659036
y(23) 	= 0.1187281829; % carel = 0.1187281829
y(24) 	= 1.5768287365; % caup = 1.5768287365
y(25) 	= 0.0000560497; % casub = 0.0000560497
y(26) 	= 0.0063427103; % Ftc = 0.0063427103
y(27) 	= 0.1296677919; % Ftmc = 0.1296677919
y(28) 	= 0.7688656371; % Ftmm = 0.7688656371
y(29) 	= 0.0242054739; % Fcms = 0.0242054739
y(30) 	= 0.0138533048; % Fcmi = 0.0138533048
y(31) 	= 0.1203184861; % Fcq = 0.1203184861
y(32) 	= 0.0000319121; % cai = 0.0000319121
y(33) 	= 0.6107148187; % q = 0.6107148187
y(34) 	= 0.7649576191; % fca = 0.7649576191
y(35) 	= 8.1179761505; % nai = 8.1179761505
y(36) 	= 139.8854603066; % ki = 139.8854603066
y(37) 	= 0.7720290515; % resting = 0.7720290515
y(38) 	= 0.0000000760; % open = 0.0000000760
y(39) 	= 0.0000000213; % inactivated = 0.0000000213
y(40)	= 0.2162168926; % resting_inactivated = 0.2162168926

timeInt = 1;
fileID = fopen ('mouseOutput.dat', 'w+');
y0 = y;
for timet=0:DELTAT:totalTime
    parvec(1) = 1;
    tspan = [0 DELTAT];
    %set numerical accuracy options for ODE solver.
    options = odeset ('RelTol', 1e-06, 'AbsTol', 1e-06);
    [t, y] = ode15s(@(t,y) f_mouseSAN (t,y,parvec), tspan, y0);
    y0 = y(end,:);
% time, voltage, carel, caup, casub, cai, fca
    if mod(timeInt, 500) == 0
        fprintf(fileID, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', timet, y(1), y(23), y(24), y(25), y(32), y(34));
   end;
   timeInt = timeInt + 1;   
end; % end of time for loop.
    
fclose(fileID);
dumi = 1;
end % end of driver function.
%
%
%
function dydt = f_mouseSAN(t,y,parvec)
%
%
R = 8.314472;
T = 310.5;  
F = 96.4845;

% ddt 				= 0.00025;
capacitance 			= 0.025;
vcell 				= 3.0;
l_cell 				= 66.3767257;
r_cell 				= 3.792956;
vrel 				= 0.0036;
vsub 				= 0.03328117;
vup  				= 0.0348;
vi 				= 1.34671883;
Mgi 				= 2.5;
nao 				= 140.0;
cao 				= 1.8;
ko 				= 5.4;
gst 				= 0.00006;
eist 				= 17.0;
gbna 				= 0.0001215;
gbca 				= 0.000015;
gbk  				= 0.0000025;
gk1 				= 0.229*0.0039228*0.9;
gks 				= 0.000299;
ecal 				= 47.0;
kmfca 				= 0.00035;
alpha_fca 			= 0.021;
all_ica_multiplier 		= 1.0;
ecat 				= 45.0;
enattxr  			= 41.5761;
multiplier2 			= 1.0;
gsus 				= 0.00039060;
inakmax_multiplier 		= 1.85;
inakmax = inakmax_multiplier*0.077;
kmnap 				= 14.0;
kmkp 				= 1.4;
K1ni 				= 395.3;
K1no 				= 1628;
K2ni 				= 2.289;
K2no 				= 561.4;
K3ni 				= 26.44;
K3no 				= 4.663;
Kci 				= 0.0207;
Kco 				= 3.663;
Kcni 				= 26.44;
Qci 				= 0.1369;
Qco 				= 0.0;
Qn 				= 0.4315;
tdifca 				= 0.04;
Prel 				= 2.5;
Krel 				= 0.0015;
nrel 				= 2.0;
Kup 				= 0.0006;
nup 				= 1.0;
Ttr 				= 40.0;
ConcTC 				= 0.031;
ConcTMC 			= 0.062;
kfTC 				= 88.8;
kfTMC 				= 237.7;
kbTC 				= 0.446;
kbTMC 				= 0.00751;
kfTMM 				= 2.277;
kbTMM 				= 0.751;
ConcCM 				= 0.045;
kfCM 				= 237.7;
kbCM 				= 0.542;
ConcCQ 				= 10.0;
kfCQ 				= 0.534;
kbCQ 				= 0.445;
koca 				= 10.0;
kom 				= 0.06;
kica 				= 0.5;
kim 				= 0.005;
eca50sr 			= 0.45;
maxsr 				= 15.0;
minsr 				= 1.0;
hsrr 				= 2.5;
pumphill 			= 2.0;
%
FRT = F/(R*T);
%
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
%
%
%
v 			= y(1);
dst 			= y(2);
fst 			= y(3);
dt 			= y(4);
ft 			= y(5);
ikr_act 		= y(6);
ikr_inact 		= y(7);
ikr_inact2 		= y(8);
iks_act 		= y(9);
fl12 			= y(10);
dl12 			= y(11);
fl13 			= y(12);
dl13 			= y(13);
r 			= y(14);
m_ttxr 			= y(15);
h_ttxr 			= y(16);
j_ttxr 			= y(17);
m_ttxs 			= y(18);
h_ttxs 			= y(19);
j_ttxs 			= y(20);
y_1_2 			= y(21);
y_4 			= y(22);
carel 			= y(23);
caup			= y(24);
casub 			= y(25);
Ftc 			= y(26);
Ftmc 			= y(27);
Ftmm 			= y(28);
Fcms 			= y(29);
Fcmi 			= y(30);
Fcq 			= y(31);
cai 			= y(32);
q 			= y(33);
fca 			= y(34);
nai 			= y(35);
ki 			= y(36);
resting 		= y(37);
open 			= y(38);
inactivated 		= y(39);
resting_inactivated 	= y(40);
%
%
%
ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);
%
% /*Ist********************************************************************/

qa = 1.0/(1.0 + exp(-(v+67.0)/5.0));
alphaqa = 1.0/(0.15*exp(-(v)/11.0)+0.2*exp(-(v)/700.0));
betaqa  =  1.0/(16.0*exp((v)/8.0)+15.0*exp((v)/50.0));
tauqa = 1.0/(alphaqa + betaqa);
alphaqi = 0.15*1.0/(3100.0*exp((v+10.0)/13.0)+700.3*exp((v+10.0)/70.0));
betaqi =  0.15*1.0/(95.7*exp(-(v+10.0)/10.0) + 50.0*exp(-(v+10.0)/700.0)) + 0.000229/(1+exp(-(v+10.0)/5.0));
qi = alphaqi/(alphaqi + betaqi);
tauqi = 1.0/(alphaqi + betaqi);
% // dst = dst + ddt*((qa-dst)/tauqa);
ddst = ((qa-dst)/tauqa);
dydt(2) = ddst;
% // fst = fst + ddt*((qi-fst)/tauqi);
dfst = ((qi-fst)/tauqi);
dydt(3) = dfst;
ist = gst*dst*fst*(v - eist);

% /* Ib ************************************************************************/

ibna = gbna*(v - ena);
ibca = gbca*(v - eca);
ibk  =  gbk*(v - ek);
ib = (ibna + ibca + ibk);
  
% /*IK1**********************************************************************/

xk1inf = 1.0/(1.0 + exp(0.070727*(v - ek)));
ik1 = gk1*xk1inf*(ko/(ko + 0.228880))*(v - ek);

% /**ICaT Cav3.1**************************************************************/

tau_dt = 1.0/(1.068*exp((v + 26.3)/30.0) + 1.068*exp(-(v + 26.3)/30.0));
dt_inf = 1.0/(1.0+exp(-(v + 26.0)/6.0));
% // dt = dt + ddt*((dt_inf - dt)/tau_dt);
dddt = ((dt_inf - dt)/tau_dt);
dydt(4) = dddt;
tau_ft = 1.0/(0.0153*exp(-(v+61.7)/83.3)+0.015*exp((v+61.7)/15.38));
ft_inf = 1.0/(1.0+exp((v + 61.7)/5.6));
% // ft = ft + ddt*((ft_inf - ft)/tau_ft);
dft = ((ft_inf - ft)/tau_ft);
dydt(5) = dft;
icat = gcat*ft*dt*(v - ecat);          

% /*Ikr********************************************************************/

ikr_act_inf = 1.0/(1.0 + exp(-(v+21.173694)/9.757086));
tau_ikr_act = 0.699821/(0.003596*exp((v)/15.339290) + 0.000177*exp(-(v)/25.868423));
% // ikr_act = ikr_act + ddt*(ikr_act_inf-ikr_act)/tau_ikr_act;     
dikr_act = (ikr_act_inf-ikr_act)/tau_ikr_act;     
dydt(6) = dikr_act;
ikr_inact_inf = 1.0/(1.0 + exp((v+20.758474-4.0)/(19.0)));
tau_ikr_inact = 0.2+0.9*1.0/(0.1*exp(v/54.645)+0.656*exp(v/106.157));
% // ikr_inact = ikr_inact + ddt*(ikr_inact_inf - ikr_inact)/tau_ikr_inact;
dikr_inact = (ikr_inact_inf - ikr_inact)/tau_ikr_inact;
dydt(7) = dikr_inact;
ikr = gkr*ikr_act*ikr_inact*(v - ek);

% /**IKs********************************************************************/

iks_act_inf = 1.0/(1.0 + exp(-(v-20.876040)/11.852723));
tau_iks_act =  1000.0/(13.097938/(1.0 + exp(-(v-48.910584)/10.630272)) + exp(-(v)/35.316539));
% // iks_act = iks_act + ddt*(iks_act_inf - iks_act)/tau_iks_act;
diks_act = (iks_act_inf - iks_act)/tau_iks_act;
dydt(9) = diks_act;
iks = gks*iks_act*iks_act*(v - eks);

% /*ICaL*******************************************************************/

if abs(v)<=0.001 
alpha_dl  = -28.39*(v+35.0)/(exp(-(v+35.0)/2.5)-1.0)+408.173;
elseif(abs(v+35.0)<=0.001)
alpha_dl  = 70.975-84.9*v/(exp(-0.208*v)-1.0);
else % if abs(v)>0.001&&abs(v+35.0)>0.001
alpha_dl  = -28.39*(v+35.0)/(exp(-(v+35.0)/2.5)-1.0)-84.9*v/(exp(-0.208*v)-1.0);
end

if abs(v-5.0)<=0.001
 beta_dl   = 28.575;
else % if abs(v-5.0)>0.001
beta_dl   = 11.43*(v-5.0)/(exp(0.4*(v-5.0))-1.0);
end
%
tau_dl  = 2000.0/(alpha_dl +beta_dl);
dl13_inf = 1.0/(1+exp(-(v+13.5)/6.0));
fl13_inf = 1.0/(1+exp((v+35.0)/7.3));
tau_fl = (7.4 + 45.77*exp(-0.5*(v+28.1)*(v+28.1)/(11*11)));
% // dl13 = dl13 + ddt*(dl13_inf - dl13)/tau_dl;
ddl13 = (dl13_inf - dl13)/tau_dl;
dydt(13) = ddl13;
% // fl13 = fl13 + ddt*(fl13_inf - fl13)/tau_fl;
dfl13 = (fl13_inf - fl13)/tau_fl;
dydt(12) = dfl13;
dl12_inf = 1.0/(1+exp(-(v+3.0)/5.0));
fl12_inf = 1.0/(1+exp((v+36.0)/4.6));
% // dl12 = dl12 + ddt*(dl12_inf - dl12)/tau_dl;
ddl12 = (dl12_inf - dl12)/tau_dl;
dydt(11) = ddl12;
% // fl12 = fl12 + ddt*(fl12_inf - fl12)/tau_fl;
dfl12 = (fl12_inf - fl12)/tau_fl;
dydt(10) = dfl12;
fca_inf = kmfca/(kmfca+casub);
taufca = fca_inf/alpha_fca;
% // fca = fca + ddt*(fca_inf - fca)/taufca;
dfca = (fca_inf - fca)/taufca;
dydt(34) = dfca;
ical12 = gcal12*fl12*dl12*fca*(v-ecal);
ical13 = gcal13*fl13*dl13*fca*(v-ecal);
   
% /**INa**********************************************************************/

fna = (9.52e-02*exp(-6.3e-2*(v+34.4))/(1+1.66*exp(-0.225*(v+63.7))))+8.69e-2; 
m3_inf_ttxr = 1.0/(1.0 + exp(-(v+45.213705)/7.219547));
h_inf_ttxr = 1.0/(1.0 + exp(-(v+62.578120 )/(-6.084036)));
m3_inf_ttxs = 1.0/(1.0 + exp(-(v+36.097331-5.0)/5.0));
h_inf_ttxs = 1.0/(1.0 + exp((v+56.0)/3.0));
m_inf_ttxr = power(m3_inf_ttxr,0.333);
m_inf_ttxs = power(m3_inf_ttxs,0.333);
tau_m = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_h = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_j = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);
% // m_ttxs = m_ttxs + ddt*(m_inf_ttxs - m_ttxs)/tau_m;
dm_ttxs = (m_inf_ttxs - m_ttxs)/tau_m;
dydt(18) = dm_ttxs;
% // h_ttxs = h_ttxs + ddt*(h_inf_ttxs - h_ttxs)/tau_h;
dh_ttxs = (h_inf_ttxs - h_ttxs)/tau_h;
dydt(19) = dh_ttxs;
% // j_ttxs = j_ttxs + ddt*(h_inf_ttxs - j_ttxs)/tau_j;
dj_ttxs = (h_inf_ttxs - j_ttxs)/tau_j;
dydt(20) = dj_ttxs;
hs = (1.0-fna)*h_ttxs+fna*j_ttxs;
tau_mr = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_hr = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_jr = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);
% // m_ttxr = m_ttxr + ddt*(m_inf_ttxr - m_ttxr)/tau_mr;
dm_ttxr = (m_inf_ttxr - m_ttxr)/tau_mr;
dydt(15) = dm_ttxr;
% // h_ttxr = h_ttxr + ddt*(h_inf_ttxr - h_ttxr)/tau_hr;
dh_ttxr = (h_inf_ttxr - h_ttxr)/tau_hr;
dydt(16) = dh_ttxr;
% // j_ttxr = j_ttxr + ddt*(h_inf_ttxr - j_ttxr)/tau_jr;
dj_ttxr = (h_inf_ttxr - j_ttxr)/tau_jr;
dydt(17) = dj_ttxr;
hsr = (1.0-fna)*h_ttxr+fna*j_ttxr;
if abs(v)>0.005
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*(F*F/(R*T))*((exp((v-ena)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*F*((exp((v-ena)*F/(R*T))-1.0));
end
%
if abs(v)>0.005
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*(F*F/(R*T))*((exp((v-enattxr)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*F*((exp((v-enattxr)*F/(R*T))-1.0));
end
%
% /**If**************************************************************************/

y_inf = 1.0/(1.0 + exp((v+vhalf_gh)/16.3));
tau_y_1_2 = 1.5049/(exp(-(v+590.3)*0.01094)+ exp((v-85.1)/17.2));
% // y_1_2 = y_1_2 + ddt*(y_inf - y_1_2)/tau_y_1_2;
dy_1_2 = (y_inf - y_1_2)/tau_y_1_2;
dydt(21) = dy_1_2;
ihk  = 0.6167*gh*y_1_2*(v - ek);
ihna = 0.3833*gh*y_1_2*(v - ena);
ih = (ihk + ihna);

% /*Ito*************************************************************************/
q_inf = 1.0/(1.0+exp((v+49.0)/13.0));
tau_q = (6.06 + 39.102/(0.57*exp(-0.08*(v+44.0))+0.065*exp(0.1*(v+45.93))))/0.67; 
% // q = q + ddt*((q_inf-q)/tau_q);
dq = ((q_inf-q)/tau_q);
dydt(33) = dq;
r_inf = 1.0/(1.0+exp(-(v-19.3)/15.0));
tau_r = (2.75+14.40516/(1.037*exp(0.09*(v+30.61))+0.369*exp(-0.12*(v+23.84))))/0.303;
% // r = r + ddt*((r_inf-r)/tau_r);
dr = ((r_inf-r)/tau_r);
dydt(14) = dr;
ito = gto*q*r*(v-ek);
 
% /*Isus***********************************************************************/
isus = gsus*r*(v-ek);
% /*Inak***********************************************************************/
inak = inakmax*(power(ko,1.2)/(power(kmkp,1.2)+power(ko,1.2)))*(power(nai,1.3)/(power(kmnap,1.3)+power(nai,1.3)))/(1.0+exp(-(v-ena+120.0)/30.0));

% /****iNaCa*******************************************************************/

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
kcasr = maxsr - (maxsr - minsr)/(1.0 + power(eca50sr/carel,hsrr));
kosrca = koca/kcasr;
kisrca = kica*kcasr;
% // resting = resting + ddt*(kim*resting_inactivated - kisrca*casub*resting - kosrca*casub*casub*resting + kom*open);
dresting = (kim*resting_inactivated - kisrca*casub*resting - kosrca*casub*casub*resting + kom*open);
dydt(37) = dresting;
% // open    = open    + ddt*(kosrca*casub*casub*resting - kom*open - kisrca*casub*open + kim*inactivated);
dopen    = (kosrca*casub*casub*resting - kom*open - kisrca*casub*open + kim*inactivated);
dydt(38) = dopen;
% // inactivated = inactivated + ddt*(kisrca*casub*open - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated);
dinactivated = (kisrca*casub*open - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated);
dydt(39) = dinactivated;
% // resting_inactivated = resting_inactivated + ddt*(kom*inactivated - kosrca*casub*casub*resting_inactivated - kim*resting_inactivated + kisrca*casub*resting);
Jrel = ks*open*(carel - casub);
dresting_inactivated = (kom*inactivated - kosrca*casub*casub*resting_inactivated - kim*resting_inactivated + kisrca*casub*resting);
Jrel = ks*open*(carel - casub);
dydt(40) = dresting_inactivated;

Jup = Pup*(power(cai/pumpkmf,pumphill) - power(caup/pumpkmr,pumphill))/(1.0 + power(cai/pumpkmf,pumphill) + power(caup/pumpkmr,pumphill));
Jtr  = (caup - carel)/Ttr;
dFtc  = kfTC*cai*(1.0-Ftc)-kbTC*Ftc;
dFtmc = kfTMC*cai*(1.0-Ftmc-Ftmm)-kbTMC*Ftmc;
dFtmm = kfTMM*Mgi*(1.0-Ftmc-Ftmm)-kbTMM*Ftmm;
dFcms = kfCM*casub*(1.0-Fcms)-kbCM*Fcms;
dFcmi = kfCM*cai*(1.0-Fcmi)-kbCM*Fcmi;
dFcq  = kfCQ*carel*(1.0-Fcq)-kbCQ*Fcq;
% // Ftc = Ftc + ddt*dFtc;
dydt(26) = dFtc;
% // Ftmc = Ftmc + ddt*dFtmc;
dydt(27) = dFtmc;
% // Ftmm = Ftmm + ddt*dFtmm;     
dydt(28) = dFtmm;
% // Fcms = Fcms + ddt*dFcms;     
dydt(29) = dFcms;
% // Fcmi = Fcmi + ddt*dFcmi;     
dydt(30) = dFcmi;
% // Fcq = Fcq + ddt*dFcq;
dydt(31) = dFcq;        
% // casub = casub + ddt*((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);
dcasub = ((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);
dydt(25) = dcasub;
% // cai = cai + ddt*((Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmi + ConcTC*dFtc + ConcTMC*dFtmc)); 
dcai = ((Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmi + ConcTC*dFtc + ConcTMC*dFtmc)); 
dydt(32) = dcai;
% // carel = carel + ddt*(Jtr - Jrel - ConcCQ*dFcq);
dcarel = (Jtr - Jrel - ConcCQ*dFcq);
dydt(23) = dcarel;
% // caup = caup + ddt*(Jup-Jtr*vrel/vup);
dcaup = (Jup-Jtr*vrel/vup);
dydt(24) = caup;
total_current = ih+ina_ttxr+ina_ttxs+ical12+ical13+iks+ikr+ik1+ist+ib+icat+inak+isus+inaca+ito;
dvdt = - total_current/capacitance;
dvdtold = dvdt;
dydt(1) = dvdt;
% // vnew = v  + ddt*dvdt;
ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);
nai_tot = ihna+ina_ttxr+ina_ttxs+3.0*inak+3.0*inaca+ist+ibna;
ki_tot = ihk+iks+ikr+ik1+ibk-2.0*inak+isus+ito;
% // nai = nai - ddt*(nai_tot)/(F*vi);
dnai = (nai_tot)/(F*vi);
dydt(35) = dnai;
% // ki = ki - ddt*(ki_tot)/(F*vi);
dki = (ki_tot)/(F*vi);
dydt(36) = dki;
%
%
%
dydt = dydt';
end % end of RHS f_mouseSAN function.
