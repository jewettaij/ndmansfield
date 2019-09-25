SHELL = /bin/sh

FAST_LOCAL_DIRECTORY = /tmp/compile_temp

INTERNAL_LIB_PATH = ../../lib

INSTALL_PATH = .

GLOBAL_LIBRARIES = \


# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cpp .o .a
.cpp.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 



LIBS = -lm



LIB_FILES = 


INCLUDES = 


OBJECT_FILES = \
globals.o \
ndmansfield.o \
hamiltonian.o \
analysis.o \
main.o


OBJECT_SRC = \
globals.h \
err_report.h \
box_dimensions.h \
random_gen.h \
box_dimensions.h \
ndmansfield.h \
hamiltonian.h \
globals.cpp \
ndmansfield.cpp \
hamiltonian.cpp \
analysis.cpp \
main.cpp


#Default target: the ndmansfield binary
#                "ndmansfield" is the lattice polymer simulator
ndmansfield: $(OBJECT_FILES) $(LIB_FILES)
	$(ANSI_CPP) $(LFLAGS) -o ndmansfield \
	$(OBJECT_FILES) \
	$(LIBS)


GENERATED_FILES = *.o *.a core ndmansfield

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" ndmansfield
	cp -f ndmansfield $(INSTALL_PATH)

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/ndmansfield

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it from Makefile.tmp to Makefile
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS) $(INCLUDES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp

# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
