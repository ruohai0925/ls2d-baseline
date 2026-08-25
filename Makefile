CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
SRC = src/io.cpp src/godunov.cpp src/linsolve.cpp src/mg.cpp src/flow.cpp src/problems_flow.cpp src/main_flow.cpp
ls2d: $(SRC) src/*.h
	$(CXX) $(CXXFLAGS) -o $@ $(SRC)
debug: $(SRC) src/*.h
	$(CXX) -O0 -g -std=c++17 -Wall -o ls2d_dbg $(SRC)
clean:
	rm -f ls2d ls2d_dbg
