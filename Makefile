CXX = c++
CXX_FLAGS = -std=c++0x -fopenmp
OPT_FLAGS = -O3 
# OPT_FLAGS += -fopenmp -ffast-math
DEBUG_FLAGS = -g -Wall -pedantic -Wno-strict-overflow
LD_FLAGS = 

ifeq ($(NOVIZ),1)
SRCS = driver.cc
else
SRCS = visual_driver.cc
LD_FLAGS += -L/usr/X11R6/lib -lm -lpthread -lX11
SRCS +=	viz/SolutionViewer.cc
endif

SRCS +=	sim/Simulator.cc
SRCS +=	cell/Cell.cc
SRCS +=	grid/Grid.cc
SRCS +=	grid/GridLevel.cc
SRCS +=	grid/BoundaryConditions.cc

%.o: %.cc
	$(CXX) $(CXX_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS) -c $< -o $@ $(LD_FLAGS)

OBJS = $(SRCS:.cc=.o)

lbmpp: $(OBJS)
	$(CXX) $(CXX_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS) $^ -o lbmpp $(LD_FLAGS)

CLEANFILES = lbmpp
CLEANFILES += $(OBJS)
CLEANFILES += *.o[0-9][0-9]*
CLEANFILES += *.e[0-9][0-9]*
CLEANFILES += *.btr

clean:
	rm -f *~ $(CLEANFILES) driver.o visual_driver.o

# eof