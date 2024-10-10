#---------------------------------------------------------------------#
# Makefile for Monte Carlo B-W code                                   #
#                                                                     #
# C. D. Woodgate, Warwick                                        2023 #
#---------------------------------------------------------------------#

# Specify a particular compiler with "make compiler=intel", etc.
# Alternatively, these options can be specified as environment variables
# e.g. "export compiler=gfortran" can be added to $HOME/.bashrc

SYSTEM=$(shell uname -s)

# Compiler specific flags
# Note: you MUST specify a compiler option. None are specified by default.

ifeq ($(strip $(compiler)),)
  MAKECMDGOALS = error
error:
	@echo 'You need to set a value for the compiler variable'
	@echo ' eg. "make compiler=gfortran"'
	@echo 'Alternatively, you can add "export compiler=gfortran" to $$HOME/.bashrc'
endif

# Intel
ifeq ($(strip $(compiler)),intel)
  FC = mpiifort
  FFLAGS = -O3 
  FFLAGS = -O0 -warn all -check all 
  FFLAGS += -module $(OBJDIR)
  CC = icc
  CFLAGS = -O3
endif

# gfortran
ifeq ($(strip $(compiler)),gfortran)
  FC = mpifort
#  FFLAGS = -O0 -Wall -Wextra -fcheck=bounds
  FFLAGS = -O3 -Wall -Wextra -fimplicit-none
  FFLAGS += -I/usr/local/include -I$(OBJDIR) -J$(OBJDIR)
  LDFLAGS=-lgcc -lopenblas
  CC=gcc -I$(INCDIR)
  CFLAGS=-O3
endif

SRCDIR=src
BINDIR=bin
DATDIR=data
OBJDIR=obj
INCDIR=include

ifeq ($(SYSTEM),Darwin)
         FFLAGS +=-I$(shell nf-config --fflags)
         LDFLAGS += $(shell nf-config --flibs) \
                    -lnetcdf -lnetcdff
else
         FFLAGS +=$(shell nf-config --fflags)
         LDFLAGS += $(shell nf-config --flibs) 
endif

# Command to use for linking and executable
LD=$(FC)
EXE=bontewarlo.run

SRCFILES=mt19937ar.c kinds.f90 shared_data.f90 io.f90 comms.f90 write_netcdf.f90 write_xyz.f90 \
	     write_diagnostics.f90 command_line.f90 c_functions.f90 display.f90 \
         energetics.f90 analytics.f90 random_site.f90 metropolis.f90 \
	     nested_sampling.f90 tmmc.f90 wang-landau.f90 energy_spectrum.f90 initialise.f90 main.f90

OBJFILES:=$(SRCFILES:.f90=.o)
OBJFILES:=$(OBJFILES:.c=.o)

VPATH = $(SRCDIR):$(SRCDIR)/core:$(OBJDIR):$(INCDIR)

alloy: $(OBJFILES)
	$(FC) $(FFLAGS) -o $(EXE) $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS)


# Purge build files and executable
clean :
	@rm -rf $(OBJDIR) $(BINDIR) $(EXE)

# Rules for building object files
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $(OBJDIR)/$@ $<

# Rules for building object files
%.o: %.c
	$(CC) $(CFLAGS) -c -o $(OBJDIR)/$@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(OBJFILES): | $(OBJDIR)

################
# Dependencies #
################
shared_data.o: kinds.o
io.o: kinds.o shared_data.o command_line.o display.o
write_netcdf.o: kinds.o shared_data.o
write_xyz.o: shared_data.o kinds.o analytics.o
write_diagnostics.o: shared_data.o kinds.o analytics.o
comms.o: kinds.o shared_data.o io.o
command_line.o: kinds.o
c_functions.o: mt19937ar.o
display.o: kinds.o shared_data.o
energetics.o: kinds.o shared_data.o c_functions.o io.o
analytics.o: shared_data.o kinds.o display.o io.o
random_site.o: shared_data.o kinds.o c_functions.o analytics.o
metropolis.o: kinds.o shared_data.o c_functions.o energetics.o random_site.o analytics.o initialise.o
nested_sampling.o: kinds.o shared_data.o c_functions.o energetics.o random_site.o analytics.o initialise.o metropolis.o
initialise.o: kinds.o shared_data.o c_functions.o energetics.o random_site.o comms.o
main.o: initialise.o shared_data.o kinds.o c_functions.o write_netcdf.o\
	write_xyz.o write_diagnostics.o command_line.o display.o metropolis.o
