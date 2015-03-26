# APLCON C++11 wrapper

This is a C++11 wrapper class around the Constrained Least Squares
fitter APLCON by Volker Blobel. Please
[see his DESY homepage how APLCON works](http://www.desy.de/~blobel/wwwcondl.html).
Also, if you use this wrapper in a publication, don't forget to cite
Blobel's original work.

Have a look at the APLCON++ example how to use this interface:

[APLCON_example.cc](src/APLCON_example.cc)

Please report bugs and improvements! The interface using variables
defined by string is considered stable, but there will be some
interface for some more complex data types available soon.

To build this project, execute
```
mkdir build
cd build
cmake ..
make
```

