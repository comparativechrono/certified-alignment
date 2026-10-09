CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra

.PHONY: all test benchmark clean
all: build/certalign

build/certalign: src/aligner.cpp include/engine.hpp include/predictor.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -Iinclude src/aligner.cpp -o $@

test: build/certalign
	python3 tests/test_correctness.py --binary build/certalign --output results/correctness.json

benchmark: build/certalign
	python3 benchmarks/run_bench.py --binary build/certalign --output-dir results/benchmark

clean:
	rm -rf build
