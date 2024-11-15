Target  = run.exe
all:$(Target)
# To include headers, use -I$PATHTOFHEADERFOLDER
# To include library, use -L$PATHTOFLIBFOLDER, and use -l$NAMEOFLIBRARY (*.so file, without the prefix lib and without the extension .so)
run.exe : tag.o
	g++ -std=c++1y tag.o -o run.exe `root-config --cflags --glibs` -lTreePlayer -lTMVA
tag.o: tag.C
	g++ -std=c++1y -o tag.o -c tag.C `root-config --cflags --glibs` -lTreePlayer
