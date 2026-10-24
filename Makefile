CXX = c++
RAYLIB = $(shell brew --prefix raylib)
DOCTEST = $(shell brew --prefix doctest)

SRC = $(wildcard src/*.cpp)

CXXFLAGS = -std=c++23 -Isrc -I$(RAYLIB)/include -g -O3 -fsanitize=address,undefined -fno-omit-frame-pointer
LDFLAGS = -L$(RAYLIB)/lib -lraylib -fsanitize=address,undefined

PROFILE_CXX = g++-16
PROFILE_CXXFLAGS = -std=c++23 -Isrc -I$(RAYLIB)/include -g -O2 -pg
PROFILE_LDFLAGS = -L$(RAYLIB)/lib -lraylib -pg

TEST_CXXFLAGS = -std=c++23 -Isrc -I$(DOCTEST)/include -g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer

# Optimized, no sanitizers: use this for timing, not for correctness.
BENCH_CXXFLAGS = -std=c++23 -Isrc -I$(DOCTEST)/include -O3 -DNDEBUG
TC ?= *timing*

sim: main.cpp $(SRC) $(wildcard src/*.hpp)
	$(CXX) $(CXXFLAGS) main.cpp $(SRC) -o sim $(LDFLAGS)

sim_profile: main.cpp $(SRC) $(wildcard src/*.hpp)
	$(PROFILE_CXX) $(PROFILE_CXXFLAGS) main.cpp $(SRC) -o sim_profile $(PROFILE_LDFLAGS)

run_tests: tests/test_main.cpp $(SRC) $(wildcard src/*.hpp)
	$(CXX) $(TEST_CXXFLAGS) tests/test_main.cpp $(SRC) -o run_tests -fsanitize=address,undefined

test: run_tests
	./run_tests

run_bench: tests/test_main.cpp $(SRC) $(wildcard src/*.hpp)
	$(CXX) $(BENCH_CXXFLAGS) tests/test_main.cpp $(SRC) -o run_bench

bench: run_bench
	./run_bench --no-skip -tc="$(TC)"

clean:
	rm -f sim sim_profile run_tests run_bench gmon.out profile.txt

.PHONY: test bench clean
