# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 

CPPDEFS = -Duse_netCDF


include ./mkmf.template.nyu


.DEFAULT:
	-echo $@ does not exist.
all: PLEV.exe
SRC =
OBJ =
clean: neat
	-rm -f .PLEV.exe.cppdefs $(OBJ) PLEV.exe
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
PLEV.exe: $(OBJ) 
	$(LD) $(OBJ) -o PLEV.exe  $(LDFLAGS)
