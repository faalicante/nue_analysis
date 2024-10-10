Target  = run.exe
all:$(Target)
# To include headers, use -I$PATHTOFHEADERFOLDER
# To include library, use -L$PATHTOFLIBFOLDER, and use -l$NAMEOFLIBRARY (*.so file, without the prefix lib and without the extension .so)
run.exe : tag_shift.o
	g++ -std=c++1y tag_shift.o -o run.exe `root-config --cflags --glibs` -lTreePlayer -lTMVA
tag_shift.o: tag_shift.C
	g++ -std=c++1y -o tag_shift.o -c tag_shift.C `root-config --cflags --glibs` -lTreePlayer
