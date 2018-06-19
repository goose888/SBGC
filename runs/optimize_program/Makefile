PROG =	FSP_PROG

SRCS =	exp_userfunc_prog.f90 main.f90 ffsqp.f qld.f

OBJS =	exp_userfunc_prog.o main.o ffsqp.o qld.o

LIBS =	

CC = icc
CFLAGS = -O0 -traceback
FC = ifort
FFLAGS = -O0 -traceback
F90 = ifort
F90FLAGS = -O0 -traceback
LDFLAGS = -s

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

main.o: exp_userfunc_prog.o
