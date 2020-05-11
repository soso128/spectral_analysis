INC=`python3 -m pybind11 --includes`

all: snrate.cpp
	g++ -O3 -Wall -shared -std=c++11 -fPIC $(INC) snrate.cpp -o snrate`python3-config --extension-suffix`

clean:
	rm -f snrate`python3-config --extension-suffix`
