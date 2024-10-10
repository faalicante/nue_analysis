Target  = run.exe
all:$(Target)
# To include headers, use -I$PATHTOFHEADERFOLDER
# To include library, use -L$PATHTOFLIBFOLDER, and use -l$NAMEOFLIBRARY (*.so file, without the prefix lib and without the extension .so)
run.exe : shift.o
	g++ -std=c++1y shift.o -o run.exe `root-config --cflags --glibs` -lTreePlayer -lTMVA
shift.o: shift.C
	g++ -std=c++1y -o shift.o -c shift.C `root-config --cflags --glibs` -lTreePlayer
