# To include headers, use -I$PATHTOFHEADERFOLDER
# To include library, use -L$PATHTOFLIBFOLDER, and use -l$NAMEOFLIBRARY (*.so file, without the prefix lib and without the extension .so)
CXXFLAGS = -O3 -march=native -funroll-loops -fomit-frame-pointer -std=c++17
LIBS = `root-config --cflags --glibs` -lTreePlayer -lTMVA -lSpectrum

Target  = run.exe
all: $(Target)

run.exe: crop_bw.o
	g++ $(CXXFLAGS) crop_bw.o -o run.exe $(LIBS)

crop_bw.o: crop_bw.C
	g++ $(CXXFLAGS) -c crop_bw.C -o crop_bw.o $(LIBS)