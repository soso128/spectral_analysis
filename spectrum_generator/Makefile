INC=`/home/giampaol/miniconda3/envs/spectral/bin/python -m pybind11 --includes`

all: snrate.cpp
	g++ -O3 -Wall -shared -std=c++11 -fPIC $(INC) snrate.cpp -o snrate`/home/giampaol/miniconda3/envs/spectral/bin/python3-config --extension-suffix`

clean:
	rm -f snrate`/home/giampaol/miniconda3/envs/spectral/bin/python3-config --extension-suffix`
