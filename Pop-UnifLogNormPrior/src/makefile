SRC = pop_deconvolution_main.c pop_format_data.c pop_mcmc.c pop_hash.c\
randgen.c pop_birthdeath.c cholesky.c normal.c
OBJ = pop_deconvolution_main.o pop_format_data.o pop_mcmc.o pop_hash.o\
randgen.o pop_birthdeath.o cholesky.o normal.o
CFLAGS = -O3
LFLAGS =
LIBS = -lm

decon : ${OBJ}
	gcc -o popdecon ${OBJ} ${LIBS}

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
