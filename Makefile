CC = gcc
RM = rm -rf
USR_DIR = /home/9yelin9/.local
OMP_DIR = /opt/mpi/gcc-4.8.5/openmpi-4.1.0
CFLAGS = -g -O2 -Wall -mcmodel=medium -I../include -I$(USR_DIR)/include -I$(OMP_DIR)/include -fopenmp
LDFLAGS = -L../lib -L$(USR_DIR)/lib -L$(OMP_DIR)/lib -fopenmp 
LINKS = -lz -lm -lopenblas
OBJS = mc.o
TARGETS = mc

.PHONY: all clean dep
.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<

all : $(TARGETS)

mc : mc.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ mc.o $(LINKS)

clean :
	$(RM) *.o
	$(RM) $(TARGET)

dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 
