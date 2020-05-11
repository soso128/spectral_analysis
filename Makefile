INC=-I/usr/include/python3.6/ 
#-I/home/elhedri/.local/include/python2.7/
#INC=`python -m pybind11 --includes` -I/disk02/usr6/elhedri/miniconda3/envs/py3/include/python3.8 -I/disk02/usr6/elhedri/miniconda3/envs/py3/include

all: likes.cc
	g++ -O3 -Wall -shared -std=c++11 -fPIC $(INC) `python3-config --includes` likes.cc -o likes`python3-config --extension-suffix`

clean:
	rm -f snrate`python3-config --extension-suffix`
