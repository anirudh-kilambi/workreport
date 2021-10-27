LIBRMN=rmn
LIBOCN=-locn
LIBS=-l$(LIBRMN) $(LIBOCN)

FFLAGS = -fp-model source -no-wrap-margin
LFLAGS =
PROGRAMS = d.cstintrp
OPTOMP=
BINDIR=${ORDENV_PLAT}/bin

# if compiling for debugging
ifeq ($(DEBUG), 1)
	# := early evaluation for self reference in the if for debug
	FFLAGS:=$(FFLAGS) -g -traceback -O0 -check
	LFLAGS:=$(LFLAGS) -g -traceback -O0 -check
endif


all: $(PROGRAMS)

%.o: %.f90
	s.f90 $(FFLAGS) -c $^ -o $@

cstintrp.o: cstintrp.F90
	s.f90 $(FFLAGS) -c $< -o $@

d.cstintrp: kdtree2.o grids.o mesh.o mesh_grid.o intrp_oce.o remap_cons_engine.o intrp_scrip.o nemo_proj.o triangle_proj.o llxy.o intrp_itfc.o cstintrp_mod.o cstintrp_py.o cstintrp.o
	s.f90 $(LFLAGS) $^ $(LIBS)  -o $@

cstintrp_py.o: cstintrp_mod.o, intrp_itfc.o, grids.o 
	s.f90 $(LFLAGS) $^ $(LIBS)  -o $@




bindist: $(PROGRAMS)
	mkdir -p $(BINDIR)
	mv $(PROGRAMS) $(BINDIR)

clean:
	/bin/rm -f *.o *mod $(PROGRAMS)

cleanbin:
	-for i in $(PROGRAMS) ; do rm -f $(BINDIR)/$$i ; done
