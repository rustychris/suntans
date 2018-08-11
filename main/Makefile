#
# Makefile for SUNTANS
#
# Oliver Fringer
# Stanford University
# fringer@stanford.edu
#
# Need to define appropriate directories in Makefile.in first!
#
# Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior University. All Rights Reserved.
#
include Makefile.in

ifneq ($(MPIHOME),)
  # For the SGI Altix
  # CC = icc
  CC = $(MPIHOME)/bin/mpicc
  CXX = $(MPIHOME)/bin/mpicxx
  MPIFILE = 
  MPIDEF = 
else
  CC = gcc
	CXX =	g++
  MPIFILE = no-mpi.c
  MPIDEF = -DNOMPI
endif
OPTFLAGS = -O2
#OPTFLAGS = -O2 -mavx -xHost -ipo -static
#OPTFLAGS=-Ofast
#OPTFLAGS=-O3 -march=native
#OPTFLAGS=-O3 -mavx
#OPTFLAGS = -O0 -g

XLDFLAGS=-lX11 -lm
XINC=/usr/include
XLIBDIR = /usr/lib64

ifneq ($(TRIANGLEHOME),)
  TRIANGLEINCLUDE = -I$(TRIANGLEHOME)
  TRIANGLELIB = $(TRIANGLEHOME)/triangle.o
  TRIANGLELIBDIR = -L$(TRIANGLEHOME)
  TRIANGLELD = 
  TRIANGLESRC = triangulate.c
else
  TRIANGLEINCLUDE =
  TRIANGLELIB =
  TRIANGLELIBDIR =
  TRIANGLELD =
  TRIANGLESRC = triangulate-notriangle.c
endif

PARMETISHOME=$(CURDIR)/ext/ParMetis-3.1.1

print-%:
	@echo '$*=$($*)'

ifneq ($(PARMETISHOME),)
  PARMETISINCLUDE = -I$(PARMETISHOME)/ParMETISLib
  PARMETISLIB = $(PARMETISHOME)/libparmetis.a $(PARMETISHOME)/libmetis.a
  PARMETISLIBDIR = -L$(PARMETISHOME)
  PARMETISLD = -lparmetis -lmetis
  PARMETISSRC = partition.c
else
  PARMETISINCLUDE = 
  PARMETISLIB =
  PARMETISLIBDIR = 
  PARMETISLD =
  PARMETISSRC = partition-noparmetis.c
endif

ifneq ($(NETCDF4HOME),)
  # optionally set, so that Makefile.in can specify exact paths as
  # needed
  NETCDFLIBDIR?=-L$(NETCDF4HOME)/lib
  NETCDFINCLUDE?=-I$(NETCDF4HOME)/include
  NETCDFLD?=-lnetcdf
  NETCDFDEF= -DUSENETCDF
  NETCDFSRC= mynetcdf.c
else
  NETCDFLIBDIR = 
  NETCDFINCLUDE = 
  NETCDFLD = 
  NETCDFDEF = -DNONETCDF
  NETCDFSRC= mynetcdf-nonetcdf.c
endif

# For the Altix
#LD = $(CC) -lmpi
LD = $(CC)
LIBS = $(PARMETISLIB) $(TRIANGLELIB) $(NETCDFLD)
LIBDIR = $(PARMETISLIBDIR) $(TRIANGLELIBDIR) $(NETCDFLIBDIR)
LDFLAGS = -lm $(LIBDIR) $(LIBS)
INCLUDES = $(PARMETISINCLUDE) $(TRIANGLEINCLUDE) $(NETCDFINCLUDE) 
DEFINES = $(MPIDEF) $(NETCDFDEF)
CFLAGS = $(OPTFLAGS) $(INCLUDES) $(DEFINES)

EXEC = sun
PEXEC = sunplot
JEXEC = sunjoin

DEPFLAGS = -Y

SRCS = 	mympi.c grid.c gridio.c report.c util.c fileio.c phys.c physio.c suntans.c initialization.c memory.c \
	turbulence.c boundaries.c check.c scalars.c tvd.c timer.c profiles.c state.c tides.c \
	sources.c diffusion.c met.c averages.c age.c merge.c sendrecv.c sediments.c \
	$(TRIANGLESRC) $(PARMETISSRC) $(MPIFILE) $(NETCDFSRC)
OBJS = $(SRCS:.c=.o)

PLOTSRCS = sunplot.c fileio.c
PLOTOBJS = $(PLOTSRCS:.c=.o)

JOINSRCS = sunjoin.c mympi.c grid.c report.c util.c fileio.c phys.c initialization.c memory.c \
	turbulence.c boundaries.c check.c scalars.c tvd.c timer.c profiles.c state.c tides.c \
	sources.c diffusion.c met.c averages no-mpi.c $(TRIANGLESRC) $(PARMETISSRC) $(NETCDFSRC)
JOINOBJS = $(JOINSRCS:.c=.o)

all:	$(EXEC)

config.vars: FORCE
	echo "MPIHOME='$(MPIHOME)'" > config.vars
	echo "PARMETISHOME='$(PARMETISHOME)'" >> config.vars
	echo "TRIANGLEHOME='$(TRIANGLEHOME)'" >> config.vars
	echo "NETCDF4HOME='$(NETCDF4HOME)'" >> config.vars
	echo "NETCDFLIBDIR='$(NETCDFLIBDIR)'" >> config.vars
	echo "NETCDFLD='$(NETCDFLD)'" >> config.vars
	echo "NETCDFINCLUDE='$(NETCDFINCLUDE)'" >> config.vars

.c.o:
	$(CC) $(CFLAGS) -c $*.c

$(EXEC): $(OBJS) 
	$(LD)  -o $@ $(OBJS) $(LDFLAGS)

$(PEXEC): $(PLOTOBJS)
	$(LD) -o sunplot $(PLOTOBJS) $(XLDFLAGS) -L$(XLIBDIR)

$(JEXEC): $(JOINOBJS)
	$(LD) -o sunjoin $(JOINOBJS) $(LDFLAGS)

ifeq ($(PARMETISHOME),$(CURDIR)/ext/ParMetis-3.1.1)
# VERSION 3 -- for now use default compiler settings
PARMETIS_URL=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/OLD/ParMetis-3.1.1.tar.gz
ext/ParMetis-3.1.1.tar.gz:
	[ -d ext ] || mkdir ext
	wget -O $@ $(PARMETIS_URL)

# --no-print-directory is an obscure workaround for an issue with calling
# sub-makes.  It solves the 'no rule to make target w' problem.
parmetis: ext/ParMetis-3.1.1/libparmetis.a

ext/ParMetis-3.1.1/libparmetis.a: ext/ParMetis-3.1.1.tar.gz
	cd ext && tar xzf $(shell basename $<)
	make --no-print-directory -C ext/ParMetis-3.1.1

clean: clean-parmetis

clean-parmetis:
	-rm -r ext/ParMetis-3.1.1

$(EXEC): parmetis

$(JEXEC): parmetis

partition.o: parmetis

endif


depend:
	makedepend $(DEPFLAGS) -- $(SRCS) $(PLOTSRCS) &> /dev/null

FORCE:

clean:
	rm -f *.o

clobber:	clean
	rm -f *~ \#*\# PI* $(EXEC) $(PEXEC) $(DEPFILE)

# DO NOT DELETE THIS LINE - Dependencies are appended after it.

mympi.o: mympi.h suntans.h fileio.h mynetcdf.h
grid.o: grid.h suntans.h fileio.h mympi.h partition.h util.h initialization.h
grid.o: memory.h triangulate.h report.h timer.h
report.o: report.h mympi.h suntans.h fileio.h grid.h
util.o: grid.h suntans.h fileio.h mympi.h util.h
fileio.o: fileio.h defaults.h suntans.h
phys.o: suntans.h phys.h grid.h fileio.h mympi.h util.h initialization.h
phys.o: memory.h turbulence.h boundaries.h check.h scalars.h timer.h
phys.o: profiles.h state.h diffusion.h sources.h met.h mynetcdf.h
suntans.o: suntans.h mympi.h fileio.h grid.h phys.h report.h
initialization.o: fileio.h suntans.h initialization.h
memory.o: memory.h
turbulence.o: phys.h suntans.h grid.h fileio.h mympi.h util.h turbulence.h
turbulence.o: boundaries.h scalars.h
boundaries.o: boundaries.h suntans.h phys.h grid.h fileio.h mympi.h mynetcdf.h
check.o: check.h grid.h suntans.h fileio.h mympi.h phys.h timer.h memory.h
scalars.o: scalars.h suntans.h grid.h fileio.h mympi.h phys.h util.h tvd.h
scalars.o: initialization.h
tvd.o: suntans.h phys.h grid.h fileio.h mympi.h tvd.h util.h
timer.o: mympi.h suntans.h fileio.h timer.h
profiles.o: util.h grid.h suntans.h fileio.h mympi.h memory.h phys.h
profiles.o: profiles.h
state.o: state.h grid.h suntans.h fileio.h mympi.h phys.h
tides.o: suntans.h mympi.h fileio.h grid.h tides.h memory.h
sources.o: phys.h suntans.h grid.h fileio.h mympi.h sources.h memory.h met.h
diffusion.o: diffusion.h grid.h suntans.h fileio.h mympi.h phys.h util.h
triangulate-notriangle.o: suntans.h mympi.h fileio.h grid.h
partition-noparmetis.o: suntans.h partition.h grid.h fileio.h mympi.h
no-mpi.o: suntans.h no-mpi.h
sunplot.o: suntans.h fileio.h
fileio.o: fileio.h defaults.h suntans.h
met.o: met.h suntans.h fileio.h memory.h grid.h phys.h util.h mynetcdf.h
mynetcdf.o: mynetcdf.h suntans.h phys.h grid.h met.h boundaries.h
averages.o: averages.h phys.h grid.h met.h
