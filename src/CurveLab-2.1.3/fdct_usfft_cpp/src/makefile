include ../../makefile.opt

LIB_SRC = 	fdct_usfft.cpp	afdct_usfft.cpp ifdct_usfft.cpp	fdct_usfft_param.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	test.cpp

MEX_SRC =	fdct_usfft_mex.cpp	afdct_usfft_mex.cpp	ifdct_usfft_mex.cpp	fdct_usfft_param_mex.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libfdct_usfft.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) libfdct_usfft.a $(LIB_OBJ)
	$(RANLIB) libfdct_usfft.a

test: libfdct_usfft.a test.o
	${CXX} -o test test.o libfdct_usfft.a ${LDFLAGS}
	./test options

matlab:	libfdct_usfft.a	${MEX_SRC}
	cp ${FFTW_DIR}/fftw/.libs/libfftw.a .
	${RANLIB} libfftw.a
	${MEX} fdct_usfft_mex.cpp ${INCLUDES} libfdct_usfft.a libfftw.a
	${MEX} afdct_usfft_mex.cpp ${INCLUDES} libfdct_usfft.a libfftw.a
	${MEX} ifdct_usfft_mex.cpp ${INCLUDES} libfdct_usfft.a libfftw.a
	${MEX} fdct_usfft_param_mex.cpp ${INCLUDES} libfdct_usfft.a libfftw.a
	rm libfftw.a
	mv *_mex.mex* ../mex

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.out libfdct_usfft.a test

tags:
	etags *hpp *cpp
