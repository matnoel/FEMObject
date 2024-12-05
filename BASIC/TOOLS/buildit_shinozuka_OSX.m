% buildit
% mex -setup C++

arch = computer('arch');

mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozuka1D.cpp
mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozuka2D.cpp
mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozuka3D.cpp
mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozukaRand1D.cpp
mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozukaRand2D.cpp
mex -I/usr/local/include -I/usr/local/include/eigen3 -L/usr/local/lib CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozukaRand3D.cpp