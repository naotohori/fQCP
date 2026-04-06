# fQCP
Quaternion Characteristic Polynomial (QCP) is a method to calculate RMSD and rotation matrix for superpositioning.

This is a Fortran version of QCP translated from C.  
The original code (QCProt 1.4; 2012, October 10) was written by Douglas L. Theobald
and Pu Liu, and distirubuted under a BSD open source license at http://theobald.brandeis.edu/qcp .

## Python Usage

This package provides a Python interface with two backends:
- **Fortran (F2PY)**: Fast compiled backend (~5-37x faster). Requires `gfortran`, `meson`, and `ninja`.
- **NumPy**: Pure Python fallback. Only requires `numpy` — no compiler needed.

The backend is selected automatically: if F2PY-compiled `.so` files are present, the Fortran backend is used; otherwise it falls back to NumPy.

```python
from fQCP import calc_rmsd, calc_rotation, superimpose
import fQCP
print(fQCP.BACKEND)  # 'fortran' or 'numpy'

# All functions take (N, 3) numpy arrays
rmsd = calc_rmsd(coords1, coords2)
rmsd, mat = calc_rotation(coords1, coords2)   # mat: 4x4 homogeneous matrix
rmsd, coords2_fit = superimpose(coords1, coords2)
```

To build the Fortran backend (optional):
```
make f2py
```

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
