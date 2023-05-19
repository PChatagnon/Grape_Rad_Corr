CC              = g++ -std=c++11
CC_OBJ_FLAGS    = -c -fPIC
CC_Shared_FLAGS = -shared -Wl,-soname,libRad_Corr_Grape.so
ROOT_CFLAGS     = $(shell ${ROOTSYS}/bin/root-config --cflags)
ROOT_LIBS       = $(shell ${ROOTSYS}/bin/root-config --libs)
libRad_Corr_Grape	= libRad_Corr_Grape

all:	    Rad_Corr_Grape.cc  RadCorr.o
	    mkdir -p lib ; rm -f lib/*.so
	    $(CC) $(CC_Shared_FLAGS) -o lib/${libRad_Corr_Grape}.so.1.0.1 RadCorr.o
	    cd lib;\
	    ln -sf ${libRad_Corr_Grape}.so.1.0.1 ${libRad_Corr_Grape}.so.1; ln -sf ${libRad_Corr_Grape}.so.1.0.1 ${libRad_Corr_Grape}.so
	    cd ../;
	    $(CC) -o Rad_Corr_Grape.exe Rad_Corr_Grape.cc -I ./include -L./lib -lRad_Corr_Grape $(ROOT_CFLAGS) $(ROOT_LIBS)
	
		
RadCorr.o: src/RadiativeCorrections.cc include/RadiativeCorrections.h
	    $(CC) $(CC_OBJ_FLAGS) src/RadiativeCorrections.cc -o $@ $(ROOT_CFLAGS) -I ./include


clean:	    
	    rm -f Rad_Corr_Grape.exe *.o lib/*.so.* lib/*.so
