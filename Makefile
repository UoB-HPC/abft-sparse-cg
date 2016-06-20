CXX      = c++
CXXFLAGS = -std=gnu++11 -I . -O3 -Wall -g
LDFLAGS  = -lm

PLATFORM = $(shell uname -s)
ARCH     = $(shell uname -p)

ifeq ($(PLATFORM), Darwin)
	LDFLAGS   = -framework OpenCL
else
	CXXFLAGS += -fopenmp
	LDFLAGS   = -lOpenCL -lm -fopenmp
endif

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
CSR/CPUContext.o: CGContext.h

CSR_OBJS += CSR/OCLContext.o
CSR/OCLContext.o: CGContext.h

ifneq (,$(findstring armv7,$(ARCH)))
  CSR_OBJS += CSR/ARM32Context.o
  CSR/ARM32Context.o: CGContext.h
endif

cg-csr: $(CSR_OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)
CSR_EXES += cg-csr



BENCHMARK_SIZE=10
benchmark: benchmark-coo benchmark-csr
benchmark-coo: cg-coo
	./run_benchmark ./cg-coo -b $(BENCHMARK_SIZE)
benchmark-csr: cg-csr
	./run_benchmark ./cg-csr -b $(BENCHMARK_SIZE)

test: test-coo test-csr
test-coo: cg-coo
	./run_tests ./cg-coo
test-csr: cg-csr
	./run_tests ./cg-csr

clean:
	rm -f cg-coo cg-csr $(COO_OBJS) $(CSR_OBJS)

.PHONY: clean test
