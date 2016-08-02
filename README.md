# The density matrix renormalization group (DMRG) code for 1D quantum systems

### Packages Required
  * The LAPACK package.

  This package is used for basic linear algebra calculations. For details see the [LAPACK modules] (http://www.netlib.org/lapack/explore-html/modules.html).

  * The ARPACK package.

  This package is used for solving large sparse matrix. 

### Running DMRG

```
make filename

./filename.out parameters
```

> 'filename' does not include the file extension '.cpp'

> 'parameters' can be omitted, then default parameters are used for the calculation

One example is given for the Kitaev Chain in the file 'kitaev.cpp'. To run the example code, replace the 'filename' by 'kitaev' in the above code.
