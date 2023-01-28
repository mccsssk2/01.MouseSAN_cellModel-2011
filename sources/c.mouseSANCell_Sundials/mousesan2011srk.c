/*
Sanjay R Kharche.
Mouse SAN cell model, Sundials version.
 my standard header and #defines, as of 22 Dec. 2013.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <cvodes/cvodes_dense.h>


#define Ith(v,i)    		NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ    */
#define IJth(A,i,j) 		DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ    */
#define RTOL  	    	RCONST(1.0e-12)   /* scalar relative tolerance            */
#define ATOL        	RCONST(1.0e-6)   /* scalar absolute tolerance components */
#define MAXSTEPS    	500000
#define ZERO        	RCONST(0.0)

#define Ith(v,i)    		NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) 		DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


#define num_beats 	50*5
#define NEQS   		40                /* number of ODEs, Bondarenko 2012 model.  */
#define NUMPARAMS 	102

#define DELTAT 0.1 // ms.

/* Type: UserData
contains problem constants.
This provides the functionality of
user defined problem constants.
for CVODE
*/
typedef struct {
realtype p[NUMPARAMS];
realtype Istim;
realtype lap_val;
} *UserData;

#include "f_san.c"

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */
int main (int argc, char *argv[])
{
  realtype t, tout;
  N_Vector y_mouseS;
  void *cvode_memS;
  int iout;
  char *str;

  FILE *output;
  int init_counter = 0;
  double i_Stim = 0.0;

/* The measurement arrays */
int NOUT, i;
double outputTime = 0.0;
double pcl;

UserData dataS;
dataS = (UserData) malloc(sizeof *dataS); // now it is created.

/* Start the automation loops here */
i_Stim 		= 0.0;
outputTime 	= 0.0;

y_mouseS 	= NULL;
cvode_memS	= NULL;

pcl = 100.0; // ms.

NOUT = (int)(((num_beats+1)*pcl)/DELTAT);

/* Create serial vector of length NEQ for I.C. and abstol */
  y_mouseS = N_VNew_Serial(NEQS);


/*************************the parameters going to cvode()************************************/

dataS->p[0] = /*           R */ 8.314472; // ) // mol^-1
dataS->p[1] = /*           T  */ 310.5; // )    
dataS->p[2] = /*           F  */ 96.4846; // )   // units: C/mol
dataS->p[3] = /*           capacitance  */ 0.025; // )  // nF
dataS->p[4] = /*           vcell  */ 3.0; // )  // pL
dataS->p[5] = /*           l_cell  */ 66.3767257; // )  // microM
dataS->p[6] = /*           r_cell  */ 3.792956; // )  // microM
dataS->p[7] = /*           vrel  */ 0.0036; // ) 
dataS->p[8] = /*           vsub  */ 0.03328117; // ) 
dataS->p[9] = /*           vup   */ 0.0348; // ) 
dataS->p[10] = /*           vi  */ 1.34671883; // ) 
dataS->p[11] = /*           Mgi  */ 2.5; // ) 	 	// Internal Mg2+ concentration
dataS->p[12] = /*           nao  */ 140.0; // )        // mM
dataS->p[13] = /*           cao  */ 1.8; // )          // mM
dataS->p[14] = /*           ko  */ 5.4; // )           // mM
dataS->p[15] = /*           gst  */ 0.00006; // )  // it can be lower. Free parameter. Maltsev has 0.000075 microS
dataS->p[16] = /*           eist  */ 17.0; // )  // mV reversal of Ist
dataS->p[17] = /*           gbna  */ 0.0001215; // )          //uS 
dataS->p[18] = /*           gbca  */ 0.000015; // )          //uS 
dataS->p[19] = /*           gbk   */ 0.0000025; // )          //uS
dataS->p[20] = /*           gk1  */ 0.000808489; // )              // microS: mangoni is 0.0009 microS
dataS->p[21] = /*           gkr  */ 0.002364; // )       // micro-S   increased it temporarily to improve AP 0.8*0.002955
dataS->p[22] = /*           gks  */ 0.000299; // )     // micro-S ki=140 mM
dataS->p[23] = /*           gcal12  */ 0.006; // ) 
dataS->p[24] = /*           gcal13  */ 0.018; // ) 
dataS->p[25] = /*           ecal  */ 47.0; // ) 
dataS->p[26] = /*           kmfca  */ 0.00035; // ) 
dataS->p[27] = /*           alpha_fca  */ 0.021; // ) 
dataS->p[28] = /*           ecat  */ 45.0; // ) 
dataS->p[29] = /*           enattxr  */ 41.5761; // ) 
dataS->p[30] = /*           gsus  */ 0.00039060; // )         // micro-S
dataS->p[31] = /*           inakmax  */ 0.14245; // )  // 2.88 pA/pF as in Kurata and Maltsev x conversion to give nA
dataS->p[32] = /*           kmnap  */ 14.0; // )  // mM AFFECTED BY ISO
dataS->p[33] = /*           kmkp  */ 1.4; // )  // mM
dataS->p[34] = /*           kNaCa  */ 5.5; // ) 
dataS->p[35] = /*           K1ni  */ 395.3; // ) 
dataS->p[36] = /*           K1no  */ 1628.0; // ) 
dataS->p[37] = /*           K2ni  */ 2.289; // ) 
dataS->p[38] = /*           K2no  */ 561.4; // ) 
dataS->p[39] = /*           K3ni  */ 26.44; // ) 
dataS->p[40] = /*           K3no  */ 4.663; // ) 
dataS->p[41] = /*           Kci  */ 0.0207; // ) 
dataS->p[42] = /*           Kco  */ 3.663; // ) 
dataS->p[43] = /*           Kcni  */ 26.44; // ) 
dataS->p[44] = /*           Qci  */ 0.1369; // ) 
dataS->p[45] = /*           Qco  */ 0.0; // ) 
dataS->p[46] = /*           Qn  */ 0.4315; // ) 
dataS->p[47] = /*           tdifca  */ 0.04; // )  // as in Maltsev. This is diffusion from Cai to Casub, this is not to do with the transients.
dataS->p[48] = /*           Prel  */ 2.5; // ) // this is useless  //	% SR Ca2+ release rate constant (ms-1): as given in the word document
dataS->p[49] = /*           Krel  */ 0.0015; // ) 
dataS->p[50] = /*           nrel  */ 2.0; // ) 	//% SR Ca2+ release Km (mM) and Hill coefficient
dataS->p[51] = /*           Kup  */ 0.0006; // )  	// 0.0006
dataS->p[52] = /*           nup  */ 1.0; // ) 		//% SR Ca2+ uptake Km (mM) and Hill coefficient
dataS->p[53] = /*           Ttr  */ 40.0; // )  	// Time constant for Ca2+ transfer from NSR to JSR (ms) 60
dataS->p[54] = /*           ConcTC  */ 0.031; // )             //              % Concentration of Troponin-Ca complex (mM)
dataS->p[55] = /*           ConcTMC  */ 0.062; // )            //      % Concentration of Troponin-Mg complex (mM)
dataS->p[56] = /*           kfTC  */ 88.8; // ) 
dataS->p[57] = /*           kfTMC  */ 237.7; // )              //% Rate constant for Ca2+ binding to Troponin(mM/ms)
dataS->p[58] = /*           kbTC  */ 0.446; // ) 
dataS->p[59] = /*           kbTMC  */ 0.00751; // )    //      % Rate constant for Ca2+ unbinding from Troponin (ms-1)
dataS->p[60] = /*           kfTMM  */ 2.277; // )              //              % Rate constant for Mg2+ binding to Troponin(mM/ms)
dataS->p[61] = /*           kbTMM  */ 0.751; // )              //              % Rate constant for Mg2+ unbinding from Troponin (ms-1)
dataS->p[62] = /*           ConcCM  */ 0.045; // )     //                      % Concentration of Calmodulin (mM)
dataS->p[63] = /*           kfCM  */ 237.7; // )               // manip.               % Rate constant for Ca2+ binding to Calmodulin (mM/ms)
dataS->p[64] = /*           kbCM  */ 0.542; // )               //              % Rate constant for Ca2+ unbinding from Calmodulin (mM/ms)
dataS->p[65] = /*           ConcCQ  */ 10.0; // )              // % Concentration of Calsequestrin (mM)
dataS->p[66] = /*           kfCQ  */ 0.534; // )               //              % Rate constant for Ca2+ binding to Calsequestrin (mM/ms)
dataS->p[67] = /*           kbCQ  */ 0.445; // )               //              % Rate constant for Ca2+ unbinding from Calsequestrin (mM/ms)
dataS->p[68] = /*           koca  */ 10.0; // )   // mM-2 ms-1
dataS->p[69] = /*           kom  */ 0.06; // )  // ms-1
dataS->p[70] = /*           kica  */ 0.5; // )  // mM-1ms-1
dataS->p[71] = /*           kim  */ 0.005; // )  // ms-1
dataS->p[72] = /*           eca50sr  */ 0.45; // )  // mM
dataS->p[73] = /*           ks  */ 1300000.0; // )  // ms-1 for a frequency of 5 hertz or more
dataS->p[74] = /*           maxsr  */ 15.0; // ) 
dataS->p[75] = /*           minsr  */ 1.0; // ) 
dataS->p[76] = /*           hsrr  */ 2.5; // ) 
dataS->p[77] = /*           pumpkmf  */ 0.000246; // )  // mM
dataS->p[78] = /*           pumpkmr  */ 3.29; // )  // mM
dataS->p[79] = /*           pumphill  */ 2.0; // )  // 2 seems to be alright
dataS->p[80] = /*           gna_ttxs  */  5.925e-06; // ) 
dataS->p[81] = /*           gna_ttxr  */  5.925e-06; // ) 
dataS->p[82] = /*           gcat  */ 0.013965; // ) 
dataS->p[83] = /*           gh  */ 0.0057; // ) 
dataS->p[84] = /*           gto ; */ 0.00492; // ) 
dataS->p[85] = /*           Pup  */ 0.04; // )
dataS->p[86] = /*              soicr_switch */ 0.0; // 1 ) // soicr usually is switched off.
dataS->p[87] = /*              soicr_max; */ 25.0; // 1 )
dataS->p[88] = /*              soicr_thresh */ 8400.0; // 1 ) // 8.4 mM, i.e. 8400 microM
dataS->p[89] = /*              soicr_slope; // 1.0; //  1 ) // 0.001 mM, i.e. 1.0 microM
dataS->p[90] = /*             soicr_tau */ 125.0; // 1 )

dataS->p[91] = /*   	   if_vhalf; */ 106.8; // ) // ISO
dataS->p[92] = /*   	   dl13vhalf */ 13.5; // )
dataS->p[93] = /*   	   fl13vhalf */ 35.0; // )
dataS->p[94] = /*   	   dl12vhalf */ 3.0; // )
dataS->p[95] = /*   	   fl12vhalf */ 36.0; 
dataS->p[96] = /*   	   ikractvhalf */ 21.173694; // )
dataS->p[97] = 1.0; //   cashiftplus;
dataS->lap_val = 0.0; //   lap_val; // this could be called lap_val, not just dataS->p[98].
dataS->p[98] = 1.0; //   cashiftminus;
dataS->p[99] = 1.0; //  ikr_act_factor;
dataS->p[100] = 0.0; //   csqnfree;
dataS->p[101] = 0.0; //  open_const;

/*************************the parameters************************************/

Ith(y_mouseS, 0+1)= -64.5216286940; // _(	v  ,  -64.5216286940)
Ith(y_mouseS, 1+1)= 0.6246780312; // _(	dst  ,  0.6246780312)
Ith(y_mouseS, 2+1)= 0.4537033169; // _(	fst  ,  0.4537033169)
Ith(y_mouseS, 3+1)= 0.0016256324; // _(	dt  ,  0.0016256324)
Ith(y_mouseS, 4+1)= 0.4264459666; // _(	ft  ,  0.4264459666)
Ith(y_mouseS, 5+1)= 0.4043600437; // _(	ikr_act  ,  0.4043600437)
Ith(y_mouseS, 6+1)= 0.9250035423; // _(	ikr_inact  ,  0.9250035423)
Ith(y_mouseS, 7+1)= 0.1875749806; // _(	ikr_inact2  ,  0.1875749806)
Ith(y_mouseS, 8+1)= 0.0127086259; // _(	iks_act  ,  0.0127086259)
Ith(y_mouseS, 9+1)= 0.9968141226; // _(	fl12  ,  0.9968141226)
Ith(y_mouseS, 10+1)= 0.0000045583; // _(	dl12  ,  0.0000045583)
Ith(y_mouseS, 11+1)= 0.9809298233; // _(	fl13  ,  0.9809298233)
Ith(y_mouseS, 12+1)= 0.0002036671; // _(	dl13  ,  0.0002036671)
Ith(y_mouseS, 13+1)= 0.0046263658; // _(	r  ,  0.0046263658)
Ith(y_mouseS, 14+1)= 0.4014088304; // _(	m_ttxr  ,  0.4014088304)
Ith(y_mouseS, 15+1)= 0.2724817537; // _(	h_ttxr  ,  0.2724817537)
Ith(y_mouseS, 16+1)= 0.0249208708; // _(	j_ttxr  ,  0.0249208708)
Ith(y_mouseS, 17+1)= 0.1079085266; // _(	m_ttxs  ,  0.1079085266)
Ith(y_mouseS, 18+1)= 0.4500098710; // _(	h_ttxs  ,  0.4500098710)
Ith(y_mouseS, 19+1)= 0.0268486392; // _(	j_ttxs  ,  0.0268486392)
Ith(y_mouseS, 20+1)= 0.0279984462; // _(	y_1_2  ,  0.0279984462)
Ith(y_mouseS, 21+1)= 0.1187281829; // _(	carel  ,  0.1187281829)
Ith(y_mouseS, 22+1)= 1.5768287365; // _(	caup  ,  1.5768287365)
Ith(y_mouseS, 23+1)= 0.0000560497; // _(	casub  ,  0.0000560497)
Ith(y_mouseS, 24+1)= 0.0063427103; // _(	Ftc  ,  0.0063427103)
Ith(y_mouseS, 25+1)= 0.1296677919; // _(	Ftmc  ,  0.1296677919)
Ith(y_mouseS, 26+1)= 0.7688656371; // _(	Ftmm  ,  0.7688656371)
Ith(y_mouseS, 27+1)= 0.0242054739; // _(	Fcms  ,  0.0242054739)
Ith(y_mouseS, 28+1)= 0.0138533048; // _(	Fcmi  ,  0.0138533048)
Ith(y_mouseS, 29+1)= 0.1203184861; // _(	Fcq  ,  0.1203184861)
Ith(y_mouseS, 30+1)= 0.0000319121; // _(	cai  ,  0.0000319121)
Ith(y_mouseS, 31+1)= 0.6107148187; // _(	q  ,  0.6107148187)
Ith(y_mouseS, 32+1)= 0.7649576191; // _(	fca  ,  0.7649576191)
Ith(y_mouseS, 33+1)= 8.1179761505; // _(	nai  ,  8.1179761505)
Ith(y_mouseS, 34+1)= 139.8854603066; // _(	ki  ,  139.8854603066)
Ith(y_mouseS, 35+1)= 0.7720290515; // _(	resting  ,  0.7720290515) // consistent initial conditions.
Ith(y_mouseS, 36+1)= 0.0000000760; // _(	open  ,  0.0000000760)
Ith(y_mouseS, 37+1)= 0.0000000213; // _(	inactivated  ,  0.0000000213)
Ith(y_mouseS, 38+1)= 0.2162168926; // _(	resting_inactivated  ,  0.2162168926)
Ith(y_mouseS, 39+1)= 0.0; // soicr

  cvode_memS = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeInit(cvode_memS, f, 0.0, y_mouseS);
  CVodeSStolerances(cvode_memS, 10e-4, 10e-6);
  CVDense(cvode_memS, NEQS);
  CVodeSetMaxStep(cvode_memS, DELTAT);

iout = 0;  tout = DELTAT;

output = fopen("mouseSAN.dat","w+");

for(iout=0;iout<NOUT;iout++){

	CVodeSetUserData(cvode_memS, dataS); 					// you must tell cvode_mem about data.
	CVode(cvode_memS, tout, y_mouseS, &t, CV_NORMAL); 	// and here is where y gets updated.

	fprintf(output, "%f %f %10.10f\n", tout, Ith(y_mouseS, 1) , Ith(y_mouseS, 31) );

      tout += DELTAT;
      outputTime = outputTime + DELTAT;
      if(outputTime >= pcl) outputTime = 0.0;
} /* end of time loop */

fclose(output);



  N_VDestroy_Serial(y_mouseS);
  CVodeFree(&cvode_memS);
  free(dataS);

	
	  return(0);
} // end of main.

