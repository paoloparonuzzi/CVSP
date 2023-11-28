EXE = KSEP

OBJS = main.o global_variables.o Compact_OK.o Graph_v4.o global_functions.o DegreeCut.o binPackingProblem.o CLIQUE.o CompactModel.o ConnectivityCut.o ComponentCut.o MixedModel.o Heuristic.o BendersCut.o MixedModelBFS.o newProblemCompCut.o

SYSTEM = x86-64_linux
LIBFORMAT = static_pic
CPLEXDIR = /opt/ibm/ILOG/CPLEX_Studio1271

CPLEXLIBDIR = $(CPLEXDIR)/cplex/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CPLEXDIR)/concert/lib/$(SYSTEM)/$(LIBFORMAT)
LP_CPLEX_INCLUDE= $(CPLEXDIR)/cplex/include
LP_CONCERT_INCLUDE = $(CPLEXDIR)/concert/include
INCDIR = -I. -I$(LP_CONCERT_INCLUDE) -I$(LP_CPLEX_INCLUDE)
ADDLIBS = -L$(CONCERTLIBDIR) -lconcert -L$(CPLEXLIBDIR) -lilocplex -lcplex -lm \
          -lpthread 
LIBS      = -static -lboost_timer -lboost_system -lboost_chrono -lboost_program_options -lboost_filesystem
ADDINCFLAGS = $(INCDIR)

SRCDIR = src
VPATH = src
OBJDIR = obj

CXX = g++
CXXFLAGS = -pipe -DIL_STD -DNDEBUG -pedantic-errors -Wno-unknown-pragmas -Wno-long-long -DBONMIN_BUILD -O3 -std=c++11 # -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion   # -Ofast -march=native -fwhole-program

INCL += $(ADDINCFLAGS)

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $(OBJDIR)/$$file`"; done; \
	$(CXX) $(CXXFLAGS) $(INCL) -o $@ $$bla $(LIBS) $(ADDLIBS)	
	
clean:
	rm -rf $(EXE) $(OBJDIR)/$(OBJS)

remake:
	make clean;
	make;
	
.cpp.o:
	@mkdir -p $(OBJDIR);\
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $(OBJDIR)/$@ `test -f '$<' || echo '$(SRCDIR)/'`$<

.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $(OBJDIR)/$@ `if test -f '$<'; \
	then $(CYGPATH_W) '$<';	else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`
