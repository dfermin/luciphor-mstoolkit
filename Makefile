MSTOOLKIT = ./mstoolkit-read-only
INCLUDES = -I$(MSTOOLKIT)/include -I/usr/include/boost
LIBPATH = -L$(MSTOOLKIT) -L/usr/lib64
LIBS = -lmstoolkit -lm -ldl -lboost_regex -lboost_filesystem -lboost_regex -lpthread -lboost_thread-mt
CPPFLAGS = -O2 -DGCC -D_FILE_OFFSET_BITS=64 

OBJS = luciphor.o PSMClass.o MSProductClass.o globals.o AscoreClass.o PepXMLClass.o statsFunctions.o nonParamDist.o FLRClass.o

CPP = g++


luciphor: $(OBJS)
	$(CPP) $(CPPFLAGS) $(INCLUDES) $(OBJS) $(LIBPATH) $(LIBS) -o luciphor
	rm -f *.o
	@echo
	@echo Done
	@echo Execute ./luciphor to get command line usage statement
	@echo

PSMClass.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c ./src/PSMClass.cpp

MSProductClass.o: 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c ./src/MSProductClass.cpp

globals.o: 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c ./src/globals.cpp

AscoreClass.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c ./src/AscoreClass.cpp

PepXMLClass.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c ./src/PepXMLClass.cpp

statsFunctions.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c ./src/statsFunctions.cpp

nonParamDist.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c ./src/nonParamDist.cpp

luciphor.o: 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c ./src/luciphor.cpp

FLRClass.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c ./src/FLRClass.cpp

clean:
	rm -f *.o
