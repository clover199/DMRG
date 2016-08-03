LAPACK_PATH = "/usr/lib"
ARPACK_PATH = "/home/y/lib/ARPACK"

CC = g++

FLAGS =

OBJS = dgemm.o dgesvd.o dsaupd.o dsyev.o \
	zgemm.o zgesvd.o zheev.o  znaupd.o \
	tensor.o qtensor.o mps.o useful.o core.o functions.o operators.o

MKL = N

ifeq ($(MKL),Y)
CC = icc
FLAGS = -mkl
OBJS = dgemm_p.o dgesvd.o dsaupd.o dsyev.o \
	zgemm_p.o zgesvd.o zheev.o  znaupd.o \
	tensor.o qtensor.o mps.o useful.o core.o functions.o operators.o
endif

MY_SRC = ./src
MY_HEAD = -I ./include/ -I ./

PATH_OBJS = $(addprefix $(MY_SRC)/, $(OBJS))

all: 
	$(MAKE) -C $(MY_SRC) CC=$(CC) FLAGS=$(FLAGS) $(OBJS)
	$(CC) $(FLAGS) test.cpp $(PATH_OBJS) $(MY_HEAD) -L$(LAPACK_PATH) -llapack -L$(ARPACK_PATH) -larpack -lblas

%: %.cpp
	$(MAKE) -C $(MY_SRC) CC=$(CC) FLAGS=$(FLAGS) $(OBJS)
	$(CC) $(FLAGS) $< $(PATH_OBJS) $(MY_HEAD) -o $@.out -L$(LAPACK_PATH) -llapack -L$(ARPACK_PATH) -larpack -lblas

clean:
	rm -rf *.o *.txt
	$(MAKE) -C ./src clean
	
fullclean:
	rm -rf *.o *.txt *.out out*
	$(MAKE) -C ./src clean