CXX = c++
CXX_FLAGS = -std=c++11
OPT_FLAGS = -O3 
# OPT_FLAGS += -ffast-math -fopenmp
DEBUG_FLAGS = -g -Wall -pedantic
LD_FLAGS = 

ifeq ($(NOVIZ),1)
SRCS = driver.cc
else
SRCS = visual_driver.cc
LD_FLAGS += -L/usr/X11R6/lib -lm -lpthread -lX11
endif

SRCS +=	src/ControlPanel.cc
SRCS +=	src/Cell.cc

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