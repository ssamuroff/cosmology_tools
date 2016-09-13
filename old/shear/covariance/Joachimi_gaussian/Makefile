include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS= -I ${GSL_INC}
USER_LDFLAGS= -lcosmosis -lm -lgslcblas 

all: interface.so

interface.so: interface.c
	$(CC) $(CFLAGS)  -shared -o interface.so interface.c $(LDFLAGS)

clean:
	rm -rf interface.so interface.so.dSYM

test:
	@echo Alas, shear/covariance has no tests

