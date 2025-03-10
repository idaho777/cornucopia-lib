This library was updated to run with Qt 6.  This has only been tested for MacOS
and not Windows.


----------
OLD README
----------
This is the source of the Cornucopia library.  This library is intended
for developers who need to turn a mouse or tablet sketch stroke into a
smooth curve.  The basic algorithm is described in:
Ilya Baran, Jaakko Lehtinen, Jovan Popovic
"Sketching Clothoid Splines Using Shortest Paths",
Eurographics 2010.

---------
LICENSING
---------

All of the source is distributed under the GNU GPL.  If you would
like to use it under a different license, contact me at
baran37@gmail.com and I will likely grant an exemption.  If you
use the library for research, please cite the above paper.

--------
BUILDING
--------

Cornucopia itself requires Eigen 3 (http://eigen.tuxfamily.org/).
As of 11/21/2010, the latest Eigen development build should work.
The demonstration UI (DemoUI) was tested with Qt 4.6.  There are
no other dependencies.

The library was tested with GCC and Visual C++ 2008 and 2010.
Both 32 and 64 bits should work.

The meta-build system is CMake, so standard instructions apply:
make a separate build directory, run cmake from it and then use
your build system.

-----
USING
-----

The interface is simple.  The example to get started is:
Test/EndToEndTest.cpp.  You control the algorithm using the
Parameters object that you pass in.  See the Parameters.h file.
If you do not want to introduce a dependency on Eigen into your
application, use the API in SimpleAPI.h.  Its use is also
demonstrated in EndToEndTest.cpp.  Oversketching is not implemented
in SimpleAPI yet--you need to call fitter.setOversketchBase,
passing the curve being oversketched.

