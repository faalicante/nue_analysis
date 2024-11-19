Target  = run.exe
all:$(Target)
# To include headers, use -I$PATHTOFHEADERFOLDER
# To include library, use -L$PATHTOFLIBFOLDER, and use -l$NAMEOFLIBRARY (*.so file, without the prefix lib and without the extension .so)
run.exe : shift_tag.o
	g++ -std=c++1y shift_tag.o -o run.exe `root-config --cflags --glibs` -lTreePlayer -lTMVA
shift_tag.o: shift_tag.C
	g++ -std=c++1y -o shift_tag.o -c shift_tag.C `root-config --cflags --glibs` -lTreePlayer
