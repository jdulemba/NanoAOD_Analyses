c++ -O3 -Wall -shared -std=c++11 -fPIC $(python -m pybind11 --includes) *.cpp -o pynusolver.so
