CC      = cc
CFLAGS  = -std=gnu99 -O3 -Wall -g
LDFLAGS = -lm

OBJS = mmio.o

ARCH = $(shell uname -p)

all: coo csr
	make -C matrices

cg-coo.o: cg.c COO/common.h
	$(CC) $(CFLAGS) -DCOO=1 -c $< -o $@
cg-csr.o: cg.c CSR/common.h
	$(CC) $(CFLAGS) -DCSR=1 -c $< -o $@


COO_OBJS = cg-coo.o COO/common.o mmio.o
COO/common.o: COO/common.h

define COO_EXE
$(1): $(2) cg-coo.o mmio.o COO/common.o
	$(CC) $$^ -o $$@ $(LDFLAGS)
$(2): COO/common.h COO/ecc.h
COO_EXES += $(1)
COO_OBJS += $(2)
endef

$(eval $(call COO_EXE, cg-coo-c, COO/spmv-c.o))
ifneq (,$(findstring armv7,$(ARCH)))
  $(eval $(call COO_EXE, cg-coo-arm32, COO/spmv-arm32.o))
endif


CSR_OBJS = cg-csr.o CSR/common.o mmio.o
CSR/common.o: CSR/common.h

define CSR_EXE
$(1): $(2) cg-csr.o mmio.o CSR/common.o
	$(CC) $$^ -o $$@ $(LDFLAGS)
$(2): CSR/common.h CSR/ecc.h
CSR_EXES += $(1)
CSR_OBJS += $(2)
endef

$(eval $(call CSR_EXE, cg-csr-c, CSR/spmv-c.o))
ifneq (,$(findstring armv7,$(ARCH)))
  $(eval $(call CSR_EXE, cg-csr-arm32, CSR/spmv-arm32.o))
endif


coo: $(COO_EXES)
csr: $(CSR_EXES)

test: coo csr
	for exe in $(COO_EXES) $(CSR_EXES); do ./run_tests ./$$exe ; done

clean:
	rm -f $(COO_EXES) $(CSR_EXES) $(COO_OBJS) $(CSR_OBJS) mmio.o

.PHONY: clean test
