CFLAGS_BASE = -O3 -mtune=native -Wno-comment -Wl,--no-as-needed -L./lib/ `root-config --cflags` `root-config --glibs` -lMinuit
CFLAGS_OPT = -g
CFLAGS_BASE += $(CFLAGS_OPT)
INCLUDE = -I./src/

CFLAGS = $(CFLAGS_BASE) $(goptical_CPPFLAGS)

LIBLOC = ./lib
#OUT = ./hdfastdirc_sim
OUT = ./hdfastdirc
#OUT = ./hdfastdirc_treeInput
#OUT = ./hdfastdirc_experiments

OBJFILES = dirc_base_sim.o
OBJFILES += dirc_threesegbox_sim.o
OBJFILES += dirc_rect_digitizer.o
OBJFILES += dirc_spread_gaussian.o
OBJFILES += GlueXUserOptions.o
#OBJFILES += DrcHit.o
#OBJFILES += DrcEvent.o

DRCLIB = /home/yunjiey/Documents/gluex_install/gluex_top/halld_recon/Linux_Ubuntu16.04-x86_64-gcc5.4.0/plugins/pid_dirc.so
#DRCLIB = /home/yunjiey/Documents/gluex_install/gluex_top/halld_recon/Linux_Ubuntu16.04-x86_64-gcc5.4.0/plugins/dirc_tree.so

#OBJFILES += dirc_optical_sim.o
#OBJFILES += dirc_babar_sim.o
#OBJFILES += dirc_lut_enum.o
#OBJFILES += dirc_lut.o
#OBJFILES += dirc_gluex_lut_enum.o
#OBJFILES += dirc_babar_digitizer.o
#OBJFILES += dirc_probability_spread.o
#OBJFILES += dirc_spread_relative.o
#OBJFILES += dirc_spread_radius.o
#OBJFILES += dirc_spread_linear_soft.o
#OBJFILES += dirc_probability_separation.o
#OBJFILES += dirc_progressive_separation.o

OBJLOC = $(patsubst %,$(LIBLOC)/%,$(OBJFILES))
LIBFILES = $(LIBLOC)
vpath %.o ./lib/
vpath %.cpp ./src/

%.o : %.cpp
	g++ -Wall $(CFLAGS) $(INCLUDE) -g -o $@ -c $<
	mv $@ $(LIBLOC)

.PHONY : all
all: $(OUT).cpp $(OBJFILES)
	g++ -Wall $(OUT).cpp $(OBJLOC) $(DRCLIB) $(CFLAGS) $(INCLUDE) -o $(OUT)

.PHONY : libs
libs: $(OBJFILES)
	echo libraries built

.PHONY : clean
clean:
	rm -f lib/*.o
	rm -f $(OUT)
	
.PHONY : cleanall
cleanall:
	rm -f lib/*
	rm -f *.gcda
	rm -f $(OUT)
	
