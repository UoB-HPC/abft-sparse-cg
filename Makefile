CC      = cc
CFLAGS  = -std=c99 -O3 -Wall -g
LDFLAGS = -lm

EXES =
COMMON_OBJS = cg.o mmio.o

all:
	make -C matrices

define COO_EXE
$(1): $(2) $(COMMON_OBJS)
	$(CC) $$^ -o $$@ $(LDFLAGS)
$(2): COO/common.h COO/ecc.h
EXES += $(1)
endef

$(eval $(call COO_EXE, cg-coo-c-baseline, COO/c/spmv-baseline.o))
$(eval $(call COO_EXE, cg-coo-c-sed, COO/c/spmv-sed.o))
$(eval $(call COO_EXE, cg-coo-c-sec7, COO/c/spmv-sec7.o))
$(eval $(call COO_EXE, cg-coo-c-sec8, COO/c/spmv-sec8.o))
$(eval $(call COO_EXE, cg-coo-c-secded, COO/c/spmv-secded.o))

all: $(EXES)

clean:
	rm -f $(EXES) $(COMMON_OBJS)

.PHONY: clean
