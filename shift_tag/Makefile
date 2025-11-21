# To include headers, use -I$PATHTOFHEADERFOLDER
# To include library, use -L$PATHTOFLIBFOLDER, and use -l$NAMEOFLIBRARY (*.so file, without the prefix lib and without the extension .so)
CXXFLAGS = -O3 -march=native -funroll-loops -fomit-frame-pointer -std=c++17
LIBS = `root-config --cflags --glibs` -lTreePlayer -lTMVA -lSpectrum

Target  = run.exe
all: $(Target)

run.exe: shift_tag.o
	g++ $(CXXFLAGS) shift_tag.o -o run.exe $(LIBS)

shift_tag.o: shift_tag.C
	g++ $(CXXFLAGS) -c shift_tag.C -o shift_tag.o $(LIBS)