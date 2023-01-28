
#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<malloc.h>
#include<time.h>

#define num_pars 25
#define if_c_or_p 1 // 1 for c and 2 for p, 3 for standard

int main(){

FILE *input, *output, *output2;
char *str;
double values[num_pars];
int i, j,imax, filecountermax,loopy, k, run;
double temp;

// printf("got here");

 if(if_c_or_p==2)
output = fopen("ptogether.dat","w");
else
output = fopen("ctogether.dat","w");

for(run=1;run<7;run++)
for(loopy = 0; loopy<1; loopy++){

  //    printf("%d %d\n",run,loopy);

str = malloc(64*sizeof(char));

if(if_c_or_p==2)
sprintf(str,"gcalGcat/Run%d/APparams.dat",run);
else
sprintf(str,"gcalGcat/Run%d/APparams.dat",run);

input = fopen(str,"r");
free(str);



for(i=0;i<num_pars;i++)
values[i] = -1.0;

 i = 0;



while(fscanf(input,"%lf",&temp)!=EOF){

values[i] = temp;

// printf("got here\n");

 if(i==7||i==10||i==13){
fprintf(output,"%10.10f\t",values[i]);
 }

i++;

 if(i==num_pars){
   if(values[13]>195&&values[13]<210&&values[24]<0.5)
     printf("%f %f %f %f %d\n",values[7],values[10],values[13],values[24],run);
   fprintf(output,"\n");
 i = 0;
 }

}

fclose(input);

}

 fclose(output);

return 0;
}
