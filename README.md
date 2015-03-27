# APLCON C++11 wrapper

This is a C++11 wrapper class around the Constrained Least Squares
fitter APLCON by Volker Blobel. Please
[see his DESY homepage how APLCON works](http://www.desy.de/~blobel/wwwcondl.html).
Also, if you use this wrapper in a publication, don't forget to cite
Blobel's original work.

Please report bugs and improvements!

## Installation and Usage

To build this project, you need
  * GNU compiler version >4.7.2 with g++ and gfortran
  * cmake version >2.6

Then do a standard out-of-source build by executing
```
mkdir build
cd build
cmake ..
make
```

You'll find some executable examples and a library called `libaplcon++.so`.

Also have a look at the provided [APLCON++ examples code](src/example)
how to use this interface.

## TODO

  * Implement linked covariances also for linked variables...
  * Actually test it as a kinfitter

