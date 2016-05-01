FETCH = wget
PLATFORM = $(shell uname -s)
ifeq ($(PLATFORM), Darwin)
	FETCH = curl -O
endif

all: shallow_water1/shallow_water1.mtx

shallow_water1/shallow_water1.mtx:
	$(FETCH) http://www.cise.ufl.edu/research/sparse/MM/MaxPlanck/shallow_water1.tar.gz
	tar xf shallow_water1.tar.gz
