LAPACK_PATH = "/usr/lib"
ARPACK_PATH = "/home/y/lib/ARPACK"

MY_SRC = "./src"

all: 
	$(MAKE) -C $(MY_SRC) all
	g++ -g test.cpp -I ./include/ -I ./\
		$(MY_SRC)/dgemm.o \
		$(MY_SRC)/dgesvd.o \
		$(MY_SRC)/dsaupd.o \
		$(MY_SRC)/dsyev.o \
		$(MY_SRC)/zgemm.o \
		$(MY_SRC)/zgesvd.o \
		$(MY_SRC)/zheev.o \
		$(MY_SRC)/znaupd.o \
		$(MY_SRC)/tensor.o \
		$(MY_SRC)/qtensor.o \
		$(MY_SRC)/functions.o \
		$(MY_SRC)/useful.o \
		$(MY_SRC)/operators.o \
		$(MY_SRC)/dcore.o \
		$(MY_SRC)/zcore.o \
		-L$(LAPACK_PATH) -llapack -L$(ARPACK_PATH) -larpack -lblas

%: %.cpp
	$(MAKE) -C $(MY_SRC) all
	g++ -g $< -I ./include/ -I./ -o $@.out\
		$(MY_SRC)/dgemm.o \
		$(MY_SRC)/dgesvd.o \
		$(MY_SRC)/dsaupd.o \
		$(MY_SRC)/dsyev.o \
		$(MY_SRC)/zgemm.o \
		$(MY_SRC)/zgesvd.o \
		$(MY_SRC)/zheev.o \
		$(MY_SRC)/znaupd.o \
		$(MY_SRC)/tensor.o \
		$(MY_SRC)/qtensor.o \
		$(MY_SRC)/functions.o \
		$(MY_SRC)/useful.o \
		$(MY_SRC)/operators.o \
		$(MY_SRC)/dcore.o \
		$(MY_SRC)/zcore.o \
		-L$(LAPACK_PATH) -llapack -L$(ARPACK_PATH) -larpack -lblas

clean:
	rm -rf *.o *.txt
	$(MAKE) -C ./src clean
	
fullclean:
	rm -rf *.o *.txt *.out
	$(MAKE) -C ./src clean