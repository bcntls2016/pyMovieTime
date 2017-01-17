COMP = ifort
CFLAGS = -c -fast -qopenmp -parallel -unroll
LDFLAGS = -qopenmp
PROGNAME2 = params
PROGNAME3 = readwf-xy

all: $(PROGNAME2) $(PROGNAME3)
	cd ..
	ln -sf src/$(PROGNAME2) ../$(PROGNAME2)
	ln -sf src/$(PROGNAME3) ../$(PROGNAME3)

#   Fortran objects
OBJS1 = $(PROGNAME1).o
OBJS2 = $(PROGNAME2).o
objs3 = readwf-xy.o
objs4 = readwf-xz.o
objss = init_deriv_parallel.o  blockds.o  blockis.o  newder.o\
			redef.o  simps.o  titols.o  interxy.o

.SUFFIXES: .f90 .f .o

$(PROGNAME1): $(OBJS1)
	$(COMP) -o $(PROGNAME1) $(OBJS1) $(LDFLAGS)
	
$(PROGNAME2): $(OBJS2)
	$(COMP) -o $(PROGNAME2) $(OBJS2) $(LDFLAGS)

$(PROGNAME3): $(objss) $(objs3)
	$(COMP)	-o $(PROGNAME3) $(objss) $(objs3) $(LDFLAGS)

$(PROGNAME4): $(objss) $(objs4)
	$(COMP)	-o $(PROGNAME4) $(objss) $(objs4) $(LDFLAGS)

.f90.o:
	$(COMP) $(CFLAGS) -o $(@) $<;

.f.o:
	$(COMP) $(CFLAGS) -o $(@) $<;

distclean:
	rm -f *.o *.bak *.lst *.mod $(PROGNAME2) $(PROGNAME3);
	rm -f ../$(PROGNAME2) ../$(PROGNAME3);

clean:
	rm -f *.o *.bak *.lst *.mod;
