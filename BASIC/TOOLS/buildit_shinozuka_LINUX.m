% buildit
% mex -setup C++

arch = computer('arch');

mex CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozuka1D.cpp
mex CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozuka2D.cpp
mex CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozuka3D.cpp
mex CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozukaRand1D.cpp
mex CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozukaRand2D.cpp
mex CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' shinozukaRand3D.cpp