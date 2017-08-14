#make file for my ecoli population model
#paths to mosek optimzer need to be adjusted to run it on different systems

INCPATHS=-I /home/warehouse/toyetunde/mosek/7/tools/platform/linux64x86/h -I.
LIBPATHS=-L /home/warehouse/toyetunde/mosek/7/tools/platform/linux64x86/bin
MOSEKLIB=-lmosek64 -pthread
CCOPT=

LDOPT=-Wl,-rpath-link,/home/warehouse/toyetunde/mosek/7/tools/platform/linux64x86/bin -Wl,-rpath,'/home/warehouse/toyetunde/mosek/7/tools/platform/linux64x86/bin' -pthread -lc -lm
CC=g++ -m64
LD=g++ -m64
CFLAG= -std=c++11

#MAINRUN=RunEcoli
#myRun: myRun.o
         
RunEcoli: RunEcoli.cc maxBiomass.tpp maxBiomass_PF.tpp
	$(CC) $(CFLAG)  -c $(INCPATHS)  $(CCOPT) -o RunEcoli.o RunEcoli.cc 
	$(LD) $(CFLAG)  $(LIBPATHS) RunEcoli.o $(MOSEKLIB) $(LDOPT) -o RunEcoli  

#all: $(MAINRUN)
	
clean:
	rm -f *.o 

	