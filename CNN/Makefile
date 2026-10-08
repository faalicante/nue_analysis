# To include headers, use -I$PATHTOFHEADERFOLDER
# To include library, use -L$PATHTOFLIBFOLDER, and use -l$NAMEOFLIBRARY (*.so file, without the prefix lib and without the extension .so)
CXXFLAGS = -O3 -march=native -funroll-loops -fomit-frame-pointer -std=c++17
LIBS = `root-config --cflags --glibs` -lTreePlayer -lTMVA -lSpectrum

Target  = data_preparation/run.exe
all: $(Target)

data_preparation/run.exe: data_preparation/crop_bw.o
	g++ $(CXXFLAGS) data_preparation/crop_bw.o -o data_preparation/run.exe $(LIBS)

data_preparation/crop_bw.o: data_preparation/crop_bw.C
	g++ $(CXXFLAGS) -c data_preparation/crop_bw.C -o data_preparation/crop_bw.o $(LIBS)