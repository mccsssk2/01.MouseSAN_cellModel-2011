#!/bin/bash
# runs.sh
# SR Kharche. 2008.
# Jonathan D. Stott <jonathan.stott@gmail.com>
# Created: Thursday, May 29, 2008 @ 11:53
# Modified: Thursday, May 29, 2008 @ 11:56

mpicc -o san mousesan.c -lm

start=1
finish=2

for (( i = $start; i <= $finish; i++  ))
do
  ar[1]=0
  ar[2]=0
  ar[3]=0
  ar[4]=0
  ar[5]=0
  ar[6]=0
  ar[7]=0
  ar[8]=0
  ar[9]=0
  ar[10]=0
  ar[11]=0
  ar[12]=0
  ar[13]=0
  ar[14]=0


    dir=dir${i}
    mkdir -p ${dir}
        cp san ${dir} 

	ar[${i}]=1

cd ${dir}
./san ${ar[1]} ${ar[2]} ${ar[3]} ${ar[4]} ${ar[5]} ${ar[6]} ${ar[7]} ${ar[8]} ${ar[9]} ${ar[10]} ${ar[11]} ${ar[12]} ${ar[13]} ${ar[14]} &
cd -
done
