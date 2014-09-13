include ../../makefile.opt

LIB_SRC = 	fdct_wrapping.cpp	ifdct_wrapping.cpp fdct_wrapping_param.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	test.cpp

MEX_SRC =	fdct_wrapping_mex.cpp	ifdct_wrapping_mex.cpp	fdct_wrapping_param_mex.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libfdct_wrapping.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) libfdct_wrapping.a $(LIB_OBJ)
	$(RANLIB) libfdct_wrapping.a

test: libfdct_wrapping.a test.o
	${CXX} -o test test.o libfdct_wrapping.a ${LDFLAGS}
	./test options

matlab:	libfdct_wrapping.a	${MEX_SRC}
	cp ${FFTW_DIR}/fftw/.libs/libfftw.a .
	${RANLIB} libfftw.a
	${MEX} fdct_wrapping_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw.a
	${MEX} ifdct_wrapping_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw.a
	${MEX} fdct_wrapping_param_mex.cpp ${INCLUDES} libfdct_wrapping.a libfftw.a
	rm libfftw.a
	mv *_mex.mex* ../mex

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.out libfdct_wrapping.a test

tags:
	etags *hpp *cpp
