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
define COO_EXE
$(1): $(2) cg-coo.o mmio.o COO/common.o
	$(CC) $$^ -o $$@ $(LDFLAGS)
$(2): COO/common.h COO/ecc.h
COO_EXES += $(1)
COO_OBJS += $(2)
endef

$(eval $(call COO_EXE, cg-coo-c-baseline, COO/c/spmv-baseline.o))
$(eval $(call COO_EXE, cg-coo-c-constraints, COO/c/spmv-constraints.o))
$(eval $(call COO_EXE, cg-coo-c-sed, COO/c/spmv-sed.o))
$(eval $(call COO_EXE, cg-coo-c-sec7, COO/c/spmv-sec7.o))
$(eval $(call COO_EXE, cg-coo-c-sec8, COO/c/spmv-sec8.o))
$(eval $(call COO_EXE, cg-coo-c-secded, COO/c/spmv-secded.o))

ifneq (,$(findstring armv7,$(ARCH)))
  $(eval $(call COO_EXE, cg-coo-arm32-sed, COO/arm32/spmv-sed.o))
endif


CSR_OBJS = cg-csr.o CSR/common.o mmio.o
define CSR_EXE
$(1): $(2) cg-csr.o mmio.o CSR/common.o
	$(CC) $$^ -o $$@ $(LDFLAGS)
$(2): CSR/common.h CSR/ecc.h
CSR_EXES += $(1)
CSR_OBJS += $(2)
endef

$(eval $(call CSR_EXE, cg-csr-c-baseline, CSR/c/spmv-baseline.o))
$(eval $(call CSR_EXE, cg-csr-c-sed, CSR/c/spmv-sed.o))
$(eval $(call CSR_EXE, cg-csr-c-sec7, CSR/c/spmv-sec7.o))
$(eval $(call CSR_EXE, cg-csr-c-sec8, CSR/c/spmv-sec8.o))
$(eval $(call CSR_EXE, cg-csr-c-secded, CSR/c/spmv-secded.o))

ifneq (,$(findstring armv7,$(ARCH)))
  $(eval $(call CSR_EXE, cg-csr-arm32-sed, CSR/arm32/spmv-sed.o))
endif


coo: $(COO_EXES)
csr: $(CSR_EXES)

test: coo csr
	for exe in $(COO_EXES) $(CSR_EXES); do \
	  ./$$exe -b 5 >/dev/null ; \
		if [ $$? -ne 0 ]; then echo "FAILED $$exe"; else echo "passed $$exe"; fi ; \
	done \

clean:
	rm -f $(COO_EXES) $(CSR_EXES) $(COO_OBJS) $(CSR_OBJS) mmio.o

.PHONY: clean test
