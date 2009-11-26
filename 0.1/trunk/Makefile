prefix = $(shell root-config --prefix)
rootconf = $(prefix)/bin/root-config 
gslconf = $(shell gsl-config --prefix)/bin/gsl-config
globesconf= $(shell globes-config --prefix)/bin/globes-config

BZDIR = /home/epp/phraar/Blitz
local_CFLAGS = -g -O4 -Wno-deprecated
INCFLAGS= -I $(shell $(rootconf) --incdir) $(shell $(globesconf) --include) -I $(BZDIR)/include 
local_LDFLAGS:=$(shell $(rootconf) --libs) $(shell $(gslconf) --libs) $(shell $(globesconf) --libs) -L $(BZDIR)/lib -lblitz
nope_prefix = /home/epp/phraar/Documents/Programs
nope_version_number = 0.1
nope_branch = trunk
nope_root = $(nope_prefix)/nope/$(nope_version_number)/$(nope_branch)

VPATH = $(nope_root)/src
BIN = ProbEngineTester 

NOPE_SRC =  prob_engine.cpp eigenfunctions.cpp ProbEngineTester.cpp
NOPE_OBJ = $(NOPE_SRC:.cpp=.o)


$(BIN) : $(NOPE_OBJ)
	gcc $(notdir $^) -o  $@  $(LDFLAGS) $(local_LDFLAGS)
%.o : %.cpp
	gcc $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
.PHONY: clean
clean:
	rm -f $(BIN) $(NOPE_OBJ)

