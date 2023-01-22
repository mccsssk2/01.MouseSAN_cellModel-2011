#!/bin/bash
#
# Kharche 2011 model make.
CC			= mpicc
CFLAGS		=
INCLUDE	= /home/kharches/software/sundials/instdir/include
MY_APP		= san
LIB			= -L/home/kharches/software/sundials/instdir/lib

san: mousesan2011srk.c
	${CC} ${CFLAGS} -I${INCLUDE} -c mousesan2011srk.c -o mousesan2011srk.o
	${CC} ${CFLAGS} mousesan2011srk.o -I${INCLUDE} -lm ${LIB} -lsundials_cvodes -lsundials_nvecserial -o ${MY_APP}

run:
	./${MY_APP}

veryclean:
	rm -r *.dat ${MY_APP} ${MY_APPI} *.o *~

