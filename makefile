CC = gcc
CFLAGS += -std=c99 -O3
# CFLAGS += -std=c99 -O3 -DGTA_BENCH
# CFLAGS += -std=c99 -g -DGTA_DEBUG
# CFLAGS += -std=c99 -g

# DEFS = -DCMEMLIMIT

GROMACS = /usr/local/gromacs
VGRO = 5

GKUT = extern/gkut
CKUT = extern/ckut

INCLUDE = include
SRC = src
TEST = test
BUILD = build
INSTALL = /usr/local/bin

LIBS = -lm

ifeq ($(VGRO),5)
INCGRO = -I$(GROMACS)/include/ \
	-I$(GROMACS)/include/gromacs/utility \
	-I$(GROMACS)/include/gromacs/fileio \
	-I$(GROMACS)/include/gromacs/commandline \
	-I$(GROMACS)/include/gromacs/legacyheaders \
	-I$(GROMACS)/include/gromacs/topology
LINKGRO = -L$(GROMACS)/lib/x86_64-linux-gnu
LIBGRO = -lgromacs
DEFS += -DGRO_V5
else
INCGRO = -I$(GROMACS)/include/gromacs
LINKGRO = -L$(GROMACS)/lib
LIBGRO = -lgmx
endif

ifneq ($(PARALLEL),0)
CFLAGS += -fopenmp
endif

MCFLAGS ='
MCFLAGS +=$(CFLAGS)
MCFLAGS +='

.PHONY: install clean

all: g_correlate test

g_correlate: g_correlate.o correlate.o gkut
	$(CC) $(CFLAGS) -o $(BUILD)/g_correlate $(BUILD)/g_correlate.o $(BUILD)/correlate.o \
	$(GKUT)/build/gkut_io.o $(GKUT)/build/gkut_log.o $(LINKGRO) $(LIBGRO) $(LIBS)

test: test_gcorr.o correlate.o gkut
	$(CC) $(CFLAGS) -o $(BUILD)/test $(BUILD)/test_gcorr.o $(BUILD)/correlate.o \
	$(GKUT)/build/gkut_io.o $(GKUT)/build/gkut_log.o $(LINKGRO) $(LIBGRO) $(LIBS)

install: g_correlate
	install $(BUILD)/g_correlate $(INSTALL)

g_correlate.o: $(SRC)/g_correlate.c $(INCLUDE)/correlate.h $(INCLUDE)/correlate.h
	$(CC) $(CFLAGS) -o $(BUILD)/g_correlate.o -c $(SRC)/g_correlate.c \
	$(DEFS) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include -I$(CKUT)/include

correlate.o: $(SRC)/correlate.c $(INCLUDE)/correlate.h
	$(CC) $(CFLAGS) -o $(BUILD)/correlate.o -c $(SRC)/correlate.c \
	$(DEFS) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include -I$(CKUT)/include

test_gcorr.o: $(TEST)/test_gcorr.c $(TEST)/minunit.h
	$(CC) $(CFLAGS) -o $(BUILD)/test_gcorr.o -c $(TEST)/test_gcorr.c \
	$(DEFS) -I$(INCLUDE) $(INCGRO) -I$(GKUT)/include -I$(CKUT)/include

gkut:
	make CC=$(CC) CFLAGS=$(MCFLAGS) GROMACS=$(GROMACS) VGRO=$(VGRO) -C $(GKUT)

clean:
	make clean -C $(GKUT) \
	&& rm -f $(BUILD)/*.o $(BUILD)/g_correlate $(BUILD)/test
