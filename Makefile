CC      = cc
CFLAGS  = -std=c99 -O3 -Wall -g
LDFLAGS = -lm

EXES = cg cg-constraints cg-sec7 cg-sec8 cg-secded
COMMON_OBJS = cg.o ecc.o mmio.o

all: $(EXES)
	cd matrices && make

cg: spmv.c common.h $(COMMON_OBJS)
	$(CC) $(CFLAGS) spmv.c      $(COMMON_OBJS) -o $@ $(LDFLAGS)

cg-constraints: spmv-constraints.c common.h $(COMMON_OBJS)
	$(CC) $(CFLAGS) spmv-constraints.c $(COMMON_OBJS) -o $@ $(LDFLAGS)

cg-sec7: spmv-sec7.c common.h $(COMMON_OBJS)
	$(CC) $(CFLAGS) spmv-sec7.c $(COMMON_OBJS) -o $@ $(LDFLAGS)

cg-sec8: spmv-sec8.c common.h $(COMMON_OBJS)
	$(CC) $(CFLAGS) spmv-sec8.c $(COMMON_OBJS) -o $@ $(LDFLAGS)

cg-secded: spmv-secded.c common.h $(COMMON_OBJS)
	$(CC) $(CFLAGS) spmv-secded.c $(COMMON_OBJS) -o $@ $(LDFLAGS)

clean:
	rm -f $(EXES) $(COMMON_OBJS)

.PHONY: clean
