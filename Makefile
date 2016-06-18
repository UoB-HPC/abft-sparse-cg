CXX      = c++
CXXFLAGS = -std=gnu++11 -I . -O3 -Wall -g
LDFLAGS  = -lm

ARCH = $(shell uname -p)

all: cg-coo cg-csr
	make -C matrices

cg-coo.o: cg.cpp COO/common.h
	$(CXX) $(CXXFLAGS) -DCOO=1 -c $< -o $@
cg-csr.o: cg.cpp CSR/common.h
	$(CXX) $(CXXFLAGS) -DCSR=1 -c $< -o $@


COO_OBJS = cg-coo.o CGContext.o mmio.o
COO_OBJS += COO/CPUContext.o
ifneq (,$(findstring armv7,$(ARCH)))
  COO_OBJS += COO/ARM32Context.o
endif

cg-coo: $(COO_OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)
COO_EXES += cg-coo


CSR_OBJS = cg-csr.o CGContext.o mmio.o
CSR_OBJS += CSR/CPUContext.o
ifneq (,$(findstring armv7,$(ARCH)))
  CSR_OBJS += CSR/ARM32Context.o
endif

cg-csr: $(CSR_OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)
CSR_EXES += cg-csr



BENCHMARK_SIZE=10
benchmark: benchmark-coo benchmark-csr
benchmark-coo:
	for exe in $(COO_EXES); do ./run_benchmark ./$$exe -b $(BENCHMARK_SIZE); done
benchmark-csr:
	for exe in $(CSR_EXES); do ./run_benchmark ./$$exe -b $(BENCHMARK_SIZE); done

test: test-coo test-csr
test-coo:
	for exe in $(COO_EXES); do ./run_tests ./$$exe ; done
test-csr:
	for exe in $(CSR_EXES); do ./run_tests ./$$exe ; done

clean:
	rm -f cg-coo cg-csr $(COO_OBJS) $(CSR_OBJS)

.PHONY: clean test
