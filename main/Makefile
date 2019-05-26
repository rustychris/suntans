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
OPTFLAGS ?= -O2

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
  substepping.c \
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
	makedepend $(DEPFLAGS) -- $(SRCS) $(PLOTSRCS) #&> /dev/null

FORCE:

clean:
	rm -f *.o

clobber:	clean
	rm -f *~ \#*\# PI* $(EXEC) $(PEXEC) $(DEPFILE)

# DO NOT DELETE THIS LINE - Dependencies are appended after it.

mympi.o: mympi.h suntans.h fileio.h mynetcdf.h phys.h grid.h met.h memory.h
mympi.o: util.h boundaries.h averages.h tvd.h age.h
grid.o: grid.h suntans.h fileio.h mympi.h partition.h util.h initialization.h
grid.o: memory.h triangulate.h report.h timer.h gridio.h sendrecv.h
gridio.o: gridio.h grid.h suntans.h fileio.h mympi.h util.h memory.h
gridio.o: sendrecv.h
report.o: report.h mympi.h suntans.h fileio.h grid.h
util.o: grid.h suntans.h fileio.h mympi.h util.h
fileio.o: fileio.h defaults.h suntans.h phys.h grid.h mympi.h
phys.o: suntans.h phys.h grid.h fileio.h mympi.h sendrecv.h util.h
phys.o: initialization.h memory.h turbulence.h boundaries.h met.h check.h
phys.o: scalars.h timer.h profiles.h state.h diffusion.h sources.h mynetcdf.h
phys.o: averages.h tvd.h age.h physio.h merge.h sediments.h substepping.h
physio.o: physio.h suntans.h phys.h grid.h fileio.h mympi.h merge.h
physio.o: sendrecv.h mynetcdf.h met.h memory.h util.h boundaries.h averages.h
physio.o: tvd.h age.h
suntans.o: suntans.h memory.h mympi.h fileio.h grid.h gridio.h phys.h
suntans.o: sediments.h boundaries.h met.h util.h physio.h report.h
initialization.o: fileio.h suntans.h initialization.h mynetcdf.h phys.h
initialization.o: grid.h mympi.h met.h memory.h util.h boundaries.h
initialization.o: averages.h tvd.h age.h
memory.o: memory.h suntans.h mympi.h fileio.h
turbulence.o: phys.h suntans.h grid.h fileio.h mympi.h sendrecv.h util.h
turbulence.o: turbulence.h boundaries.h met.h memory.h scalars.h
boundaries.o: sendrecv.h grid.h suntans.h fileio.h mympi.h boundaries.h
boundaries.o: phys.h met.h memory.h util.h mynetcdf.h averages.h tvd.h age.h
boundaries.o: sediments.h
check.o: check.h grid.h suntans.h fileio.h mympi.h phys.h timer.h memory.h
scalars.o: scalars.h suntans.h grid.h fileio.h mympi.h phys.h util.h tvd.h
scalars.o: initialization.h
tvd.o: suntans.h phys.h grid.h fileio.h mympi.h tvd.h util.h sendrecv.h
timer.o: mympi.h suntans.h fileio.h timer.h
profiles.o: util.h grid.h suntans.h fileio.h mympi.h memory.h phys.h
profiles.o: sediments.h boundaries.h met.h profiles.h
state.o: state.h grid.h suntans.h fileio.h mympi.h phys.h
tides.o: suntans.h mympi.h fileio.h grid.h tides.h memory.h
sources.o: phys.h suntans.h grid.h fileio.h mympi.h sources.h met.h memory.h
sources.o: util.h boundaries.h sendrecv.h
diffusion.o: diffusion.h grid.h suntans.h fileio.h mympi.h phys.h util.h
met.o: met.h suntans.h fileio.h memory.h grid.h mympi.h phys.h util.h
met.o: mynetcdf.h boundaries.h averages.h tvd.h age.h sendrecv.h
averages.o: averages.h suntans.h phys.h grid.h fileio.h mympi.h met.h
averages.o: memory.h util.h tvd.h age.h sendrecv.h
age.o: age.h phys.h suntans.h grid.h fileio.h mympi.h memory.h sendrecv.h
age.o: scalars.h mynetcdf.h met.h util.h boundaries.h averages.h tvd.h
merge.o: merge.h suntans.h grid.h fileio.h mympi.h memory.h
sendrecv.o: sendrecv.h grid.h suntans.h fileio.h mympi.h timer.h memory.h
sendrecv.o: util.h
sediments.o: suntans.h grid.h fileio.h mympi.h phys.h sendrecv.h timer.h
sediments.o: initialization.h boundaries.h met.h memory.h util.h sediments.h
sediments.o: mynetcdf.h averages.h tvd.h age.h scalars.h physio.h
substepping.o: suntans.h phys.h grid.h fileio.h mympi.h sendrecv.h util.h
substepping.o: timer.h substepping.h
triangulate.o: suntans.h mympi.h fileio.h grid.h triangulate.h
partition.o: partition.h grid.h suntans.h fileio.h mympi.h memory.h
mynetcdf.o: mynetcdf.h suntans.h phys.h grid.h fileio.h mympi.h met.h
mynetcdf.o: memory.h util.h boundaries.h averages.h tvd.h age.h merge.h
mynetcdf.o: sediments.h
sunplot.o: suntans.h fileio.h cmap.h
fileio.o: fileio.h defaults.h suntans.h phys.h grid.h mympi.h
