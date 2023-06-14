SRCS := $(wildcard *.f)
OBJS = $(SRCS:.f=.o)
INCS = $(wildcard *.inc)

# For Pentium IV
#FC=ifort
#FFLAGS=-xP -tpp7 -ftz -fpe0 -O3 -ip -parallel
#LIBS=-lmkl_solver -lmkl_intel -lmkl_intel_thread -lmkl_core -lguide -lpthread -lm -lvffpack

# Intel Compiler for any x86 machine
# FC=ifort
# FFLAGS=-ftz -fpe0 -O3 -ip -parallel

# Debugging for Pentium IV
# FC=ifort
# FFLAGS=-arch pn4 -xN -tpp7 -tune pn4 -ftz -fpe0 -g

# g77 for Pentium IV
FC=gfortran
FFLAGS=
LIBS= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl


# g77 for AMD Athlon MP
# FC=g77
# FFLAGS=-march=athlon-mp -O3 -mfpmath=sse -fomit-frame-pointer -ffast-math

# g77 debugging for AMD Athlong MP
# FC=g77
# FFLAGS=-march=athlon-mp -g

iim2fluid: $(SRCS) $(INCS) makefile
	$(FC) $(FFLAGS) $(SRCS) $(MKL)$ $(LIBS) -o $@
	strip $@

clean:
	rm -rf *.o


