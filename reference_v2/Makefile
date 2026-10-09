CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra

.PHONY: all test test-bounds benchmark benchmark-bounds clean
all: build/certalign

build/certalign: src/aligner.cpp include/engine.hpp include/predictor.hpp include/lower_bounds.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -Iinclude src/aligner.cpp -o $@

test: build/certalign
	python3 tests/test_correctness.py --binary build/certalign --output results/correctness.json

build/test_lower_bounds: tests/test_lower_bounds.cpp include/lower_bounds.hpp include/engine.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -Iinclude tests/test_lower_bounds.cpp -o $@

reference_v1/build/certalign: reference_v1/src/aligner.cpp reference_v1/include/engine.hpp reference_v1/include/predictor.hpp
	$(MAKE) -C reference_v1

test-bounds: build/certalign build/test_lower_bounds reference_v1/build/certalign
	mkdir -p results
	./build/test_lower_bounds > results/lower_bounds.json
	python3 tests/test_bounds_integration.py --binary build/certalign --baseline reference_v1/build/certalign --output results/bounds_integration.json

benchmark: build/certalign
	python3 benchmarks/run_bench.py --binary build/certalign --output-dir results/benchmark

benchmark-bounds: build/certalign
	python3 benchmarks/run_bound_bench.py --binary build/certalign --output-dir results/bound_benchmark

clean:
	rm -rf build
