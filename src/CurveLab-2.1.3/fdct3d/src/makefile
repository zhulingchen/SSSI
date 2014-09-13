include ../../makefile.opt

LIB_SRC = 	fdct3d_forward.cpp	fdct3d_inverse.cpp	fdct3d_param.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	test.cpp

MEX_SRC = 	fdct3d_forward_mex.cpp fdct3d_inverse_mex.cpp fdct3d_param_mex.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libfdct3d.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

test: libfdct3d.a test.o
	${CXX} -o test test.o libfdct3d.a ${LDFLAGS}
	./test options

matlab: libfdct3d.a	${MEX_SRC}
	cp ${FFTW_DIR}/fftw/.libs/libfftw.a .
	${RANLIB} libfftw.a
	${MEX} fdct3d_forward_mex.cpp ${INCLUDES} libfdct3d.a libfftw.a
	${MEX} fdct3d_inverse_mex.cpp ${INCLUDES} libfdct3d.a libfftw.a
	${MEX} fdct3d_param_mex.cpp ${INCLUDES} libfdct3d.a libfftw.a
	rm libfftw.a
	mv *_mex.mex* ../mex

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.out libfdct3d.a test

tags:
	etags *.hpp *.cpp
