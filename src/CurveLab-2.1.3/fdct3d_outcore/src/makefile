include ../../makefile.opt

LIB_SRC = 	fdct3d_param.cpp	fdct3d_forward.cpp	fdct3d_inverse.cpp	cpxcrvletocr.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	test.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libfdct3d.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

test: libfdct3d.a test.o
	${CXX} -o test test.o libfdct3d.a ${LDFLAGS}
	./test options

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.out libfdct3d.a test

tags:
	etags *.hpp *.cpp
