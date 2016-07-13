# fQCP
Quaternion Characteristic Polynomial (QCP) is a method to calculate RMSD and rotation matrix for superpositioning.

This is a Fortran version of QCP translated from C.  
The original code (QCProt 1.4; 2012, October 10) was written by Douglas L. Theobald
and Pu Liu, and distirubuted under a BSD open source license at http://theobald.brandeis.edu/qcp .

In this package, a Python module is also provided through f2py.

Following is the original notice for the original code written in C
by Theobald and Liu, distributed at http://theobald.brandeis.edu/qcp 

>    If you use this QCP rotation calculation method in a publication, please
>    reference:
>
>      Douglas L. Theobald (2005)
>      "Rapid calculation of RMSD using a quaternion-based characteristic
>      polynomial."
>      Acta Crystallographica A 61(4):478-480.
>
>      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2010)
>      "Fast determination of the optimal rotational matrix for weighted 
>      superpositions."
>      Journal of Computational Chemistry 31(7):1561-1563
>
>
>  Copyright (c) 2009-2012, Pu Liu and Douglas L. Theobald
>  All rights reserved.


## License and Copyright
Both this package and the original code are distributed under a BSD license.  
See https://github.com/naotohori/fQCP/blob/master/LICENSE

Copyright (c) 2009-2012, Pu Liu and Douglas L. Theobald  (The original code in c)  
Copyright (c) 2016, Naoto Hori  (This package)
