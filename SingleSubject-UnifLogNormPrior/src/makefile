SRC = deconvolution_main.c format_data.c mcmc.c hash.c\
randgen.c birthdeath.c cholesky.c normal.c
OBJ = deconvolution_main.o format_data.o mcmc.o hash.o\
randgen.o birthdeath.o cholesky.o normal.o
CFLAGS = -O3
LFLAGS =
LIBS = -lm

decon : ${OBJ}
	gcc -o decon ${OBJ} ${LIBS}

${OBJ} : ${SRC}
	gcc -c ${CFLAGS}  ${SRC}

#SRC = predict.c randgen.c hash.c birthdeath.c cholesky.c
#OBJ = predict.o randgen.o hash.o birthdeath.o cholesky.o

#predict : ${OBJ}
#	gcc -o predict ${OBJ} ${LIBS}

#${OBJ} : ${SRC}
#	gcc -c -O4 ${SRC}

#SRC = events.c
#OBJ = events.o

#events : ${OBJ}
#	gcc -o events ${OBJ} ${LIBS}

#${OBJ} : ${SRC}
#	gcc -c -O4 ${SRC}

{OBJ} : deconvolution_main.h

clean:
	rm ${OBJ}
