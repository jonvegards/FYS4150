# General makefile for c - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CC= c++ -Wall
PROG= myprogram
# Here we make the executable file
${PROG} :		${PROG}.o
				${CC} ${PROG}.o -o ${PROG}
# whereas here we create the object file
${PROG}.o :		${PROG}.cpp
				${CC} -c ${PROG}.cpp