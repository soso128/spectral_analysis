INC=-I/usr/include/python3.6/ 
#-I/home/elhedri/.local/include/python2.7/
#INC=`python -m pybind11 --includes` -I/disk02/usr6/elhedri/miniconda3/envs/py3/include/python3.8 -I/disk02/usr6/elhedri/miniconda3/envs/py3/include

all: likes.cc
	g++ -O3 -Wall -shared -std=c++11 -fPIC $(INC) `/home/giampaol/miniconda3/envs/spectral/bin/python -m pybind11 --includes` likes.cc -o likes`/home/giampaol/miniconda3/envs/spectral/bin/python3-config --extension-suffix`

clean:
	rm -f snrate`/home/giampaol/miniconda3/envs/spectral/bin/python3-config --extension-suffix`
