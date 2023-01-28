
#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<malloc.h>
#include<time.h>

#define num_pars 8

int main(){

FILE *input, *output, *output2;
char *str;
double values[num_pars];
int i, j,imax, filecountermax,loopy, k, run;
double temp;

output = fopen("caSubsystemData.dat","w");

 for(run=1;run<8;run++){

str = malloc(64*sizeof(char));

sprintf(str,"dir%d/Params.dat",run);

input = fopen(str,"r");
free(str);

for(i=0;i<num_pars;i++)
values[i] = -1.0;

i = 0;

while(fscanf(input,"%lf",&temp)!=EOF){

values[i] = temp;


 if(i==2||i==3||i==5){
fprintf(output,"%10.10f\t",values[i]);
 }

i++;

 if(i==num_pars){
   fprintf(output,"\n");
 i = 0;
 }

}

fclose(input);

}

 fclose(output);

return 0;
}
