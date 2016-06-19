CXX      = c++
CXXFLAGS = -std=gnu++11 -I . -O3 -Wall -g
LDFLAGS  = -lm

ARCH = $(shell uname -p)

all: cg-coo cg-csr
	make -C matrices

cg.o: CGContext.h
CGContext.o: CGContext.h


COO_OBJS = cg.o CGContext.o mmio.o

COO_OBJS += COO/CPUContext.o
COO/CPUContext.o: CGContext.h

ifneq (,$(findstring armv7,$(ARCH)))
  COO_OBJS += COO/ARM32Context.o
  COO/ARM32Context.o: CGContext.h
endif

cg-coo: $(COO_OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)
COO_EXES += cg-coo


CSR_OBJS = cg.o CGContext.o mmio.o

CSR_OBJS += CSR/CPUContext.o
COO/CPUContext.o: CGContext.h

ifneq (,$(findstring armv7,$(ARCH)))
  CSR_OBJS += CSR/ARM32Context.o
  CSR/ARM32Context.o: CGContext.h
endif

cg-csr: $(CSR_OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)
CSR_EXES += cg-csr



BENCHMARK_SIZE=10
benchmark: benchmark-coo benchmark-csr
benchmark-coo:
	./run_benchmark ./cg-coo -b $(BENCHMARK_SIZE)
benchmark-csr:
	./run_benchmark ./cg-coo -b $(BENCHMARK_SIZE)

test: test-coo test-csr
test-coo:
	./run_tests ./cg-coo
test-csr:
	./run_tests ./cg-csr

clean:
	rm -f cg-coo cg-csr $(COO_OBJS) $(CSR_OBJS)

.PHONY: clean test
