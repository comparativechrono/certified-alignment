CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra
REPEAT_OUTPUT ?= results/repeat_reproduction

.PHONY: all test test-bounds test-ordered test-sparse benchmark benchmark-bounds benchmark-ordered benchmark-sparse benchmark-repeat clean
all: build/certalign

build/certalign: src/aligner.cpp include/engine.hpp include/predictor.hpp include/lower_bounds.hpp include/ordered_bounds.hpp include/sparse_bounds.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -Iinclude src/aligner.cpp -o $@

test: build/certalign
	python3 tests/test_correctness.py --binary build/certalign --output results/correctness.json

build/test_lower_bounds: tests/test_lower_bounds.cpp include/lower_bounds.hpp include/engine.hpp include/ordered_bounds.hpp include/sparse_bounds.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -Iinclude tests/test_lower_bounds.cpp -o $@

reference_v1/build/certalign: reference_v1/src/aligner.cpp reference_v1/include/engine.hpp reference_v1/include/predictor.hpp
	$(MAKE) -C reference_v1

test-bounds: build/certalign build/test_lower_bounds reference_v1/build/certalign
	mkdir -p results
	./build/test_lower_bounds > results/lower_bounds.json
	python3 tests/test_bounds_integration.py --binary build/certalign --baseline reference_v1/build/certalign --output results/bounds_integration.json

build/test_ordered_bounds: tests/test_ordered_bounds.cpp include/lower_bounds.hpp include/ordered_bounds.hpp include/engine.hpp include/sparse_bounds.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -Iinclude tests/test_ordered_bounds.cpp -o $@

reference_v2/build/certalign: reference_v2/src/aligner.cpp reference_v2/include/engine.hpp reference_v2/include/predictor.hpp reference_v2/include/lower_bounds.hpp
	$(MAKE) -C reference_v2

test-ordered: build/certalign build/test_ordered_bounds reference_v2/build/certalign
	mkdir -p results
	./build/test_ordered_bounds > results/ordered_bounds.json
	python3 tests/test_ordered_integration.py --binary build/certalign --baseline reference_v2/build/certalign --output results/ordered_integration.json

build/test_sparse_bounds: tests/test_sparse_bounds.cpp include/sparse_bounds.hpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -Iinclude tests/test_sparse_bounds.cpp -o $@

reference_v3/build/certalign: reference_v3/src/aligner.cpp reference_v3/include/engine.hpp reference_v3/include/predictor.hpp reference_v3/include/lower_bounds.hpp reference_v3/include/ordered_bounds.hpp
	$(MAKE) -C reference_v3

test-sparse: build/certalign build/test_sparse_bounds reference_v3/build/certalign
	mkdir -p results
	./build/test_sparse_bounds > results/sparse_bounds.json
	python3 tests/test_sparse_integration.py --binary build/certalign --baseline reference_v3/build/certalign --output results/sparse_integration.json

benchmark: build/certalign
	python3 benchmarks/run_bench.py --binary build/certalign --output-dir results/benchmark

benchmark-bounds: build/certalign
	python3 benchmarks/run_bound_bench.py --binary build/certalign --output-dir results/bound_benchmark

benchmark-ordered: build/certalign
	python3 benchmarks/run_ordered_bench.py --binary build/certalign --output-dir results/ordered_benchmark

benchmark-sparse: build/certalign
	python3 benchmarks/run_sparse_bench.py --binary build/certalign --output-dir results/sparse_benchmark

build/reference_distance: tests/reference_distance.cpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) tests/reference_distance.cpp -o $@

benchmark-repeat: build/certalign build/reference_distance
	python3 benchmarks/run_repeat_bench.py --binary build/certalign --output-dir $(REPEAT_OUTPUT) --allow-rebuilt-binary --generate-only
	python3 tests/validate_repeat_oracle.py --binary build/reference_distance --pairs $(REPEAT_OUTPUT)/pairs.tsv --metadata $(REPEAT_OUTPUT)/cases.json --output $(REPEAT_OUTPUT)/oracle_scores.json
	python3 benchmarks/run_repeat_bench.py --binary build/certalign --output-dir $(REPEAT_OUTPUT) --allow-rebuilt-binary

clean:
	rm -rf build
