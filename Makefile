CC      = cc
CFLAGS  = -std=c99 -O3 -Wall -g
LDFLAGS = -lm

EXES = cg-baseline cg-constraints cg-sed cg-sec7 cg-sec8 cg-secded
COMMON_OBJS = cg.o ecc.o mmio.o

all: $(EXES)
	cd matrices && make

cg%: spmv%.o common.h $(COMMON_OBJS)
	$(CC) $< $(COMMON_OBJS) -o $@ $(LDFLAGS)

clean:
	rm -f $(EXES) $(COMMON_OBJS)

.PHONY: clean
