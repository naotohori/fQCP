"""
Pure Python/NumPy implementation of the QCP (Quaternion Characteristic Polynomial)
algorithm for rapid RMSD calculation and optimal superposition.

Based on the Fortran implementation by Naoto Hori, which is a translation of the
original C code by Douglas L. Theobald and Pu Liu.

References:
    Theobald, D. L. (2005) Acta Crystallographica A 61(4):478-480.
    Liu, P., Agrafiotis, D. K., & Theobald, D. L. (2009) J. Comput. Chem. 31(7):1561-1563.
"""

import numpy as np

_EVALPREC = 1.0e-11
_EVECPREC = 1.0e-6


def _center_coords(coords):
    """Center coords (N,3) and return (centered_copy, centroid)."""
    centroid = coords.mean(axis=0)
    return coords - centroid, centroid


def _inner_product(c1, c2):
    """Compute inner product matrix A(9) and E0 from centered (N,3) arrays."""
    G = np.sum(c1 * c1) + np.sum(c2 * c2)
    E0 = G * 0.5
    # A[0:3] = sum over i of c1[i,0]*c2[i,:] etc.
    A = np.empty(9)
    A[0] = np.dot(c1[:, 0], c2[:, 0])
    A[1] = np.dot(c1[:, 0], c2[:, 1])
    A[2] = np.dot(c1[:, 0], c2[:, 2])
    A[3] = np.dot(c1[:, 1], c2[:, 0])
    A[4] = np.dot(c1[:, 1], c2[:, 1])
    A[5] = np.dot(c1[:, 1], c2[:, 2])
    A[6] = np.dot(c1[:, 2], c2[:, 0])
    A[7] = np.dot(c1[:, 2], c2[:, 1])
    A[8] = np.dot(c1[:, 2], c2[:, 2])
    return A, E0


def _calc_characteristic_polynomial(A):
    """Compute the three coefficients C[0..2] of the characteristic polynomial."""
    Sxx, Sxy, Sxz = A[0], A[1], A[2]
    Syx, Syy, Syz = A[3], A[4], A[5]
    Szx, Szy, Szz = A[6], A[7], A[8]

    Sxx2 = Sxx * Sxx; Syy2 = Syy * Syy; Szz2 = Szz * Szz
    Sxy2 = Sxy * Sxy; Syz2 = Syz * Syz; Sxz2 = Sxz * Sxz
    Syx2 = Syx * Syx; Szy2 = Szy * Szy; Szx2 = Szx * Szx

    SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    C2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C1 = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx
                - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz)

    SxzpSzx = Sxz + Szx; SyzpSzy = Syz + Szy; SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy; SxzmSzx = Sxz - Szx; SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy; SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

    C0 = (Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz)))

    return C0, C1, C2


def _newton_raphson_max_eigenvalue(C0, C1, C2, E0):
    """Find the maximum eigenvalue via Newton-Raphson iteration."""
    mxEigenV = E0
    for _ in range(50):
        oldg = mxEigenV
        x2 = mxEigenV * mxEigenV
        b = (x2 + C2) * mxEigenV
        a = b + C1
        delta = (a * mxEigenV + C0) / (2.0 * x2 * mxEigenV + b + a)
        mxEigenV -= delta
        if abs(mxEigenV - oldg) < abs(_EVALPREC * mxEigenV):
            return mxEigenV
    return mxEigenV


def _fast_calc_rmsd(A, nlen, E0):
    """Calculate RMSD only (no rotation matrix)."""
    C0, C1, C2 = _calc_characteristic_polynomial(A)
    mxEigenV = _newton_raphson_max_eigenvalue(C0, C1, C2, E0)
    return np.sqrt(abs(2.0 * (E0 - mxEigenV) / nlen))


def _fast_calc_rmsd_and_rotation(A, nlen, E0):
    """Calculate RMSD and 3x3 rotation matrix (flattened to 9 elements)."""
    Sxx, Sxy, Sxz = A[0], A[1], A[2]
    Syx, Syy, Syz = A[3], A[4], A[5]
    Szx, Szy, Szz = A[6], A[7], A[8]

    C0, C1, C2 = _calc_characteristic_polynomial(A)
    mxEigenV = _newton_raphson_max_eigenvalue(C0, C1, C2, E0)
    rmsd = np.sqrt(abs(2.0 * (E0 - mxEigenV) / nlen))

    SyzmSzy = Syz - Szy; SxzmSzx = Sxz - Szx; SxymSyx = Sxy - Syx
    SxzpSzx = Sxz + Szx; SyzpSzy = Syz + Szy; SxypSyx = Sxy + Syx
    SxxpSyy = Sxx + Syy; SxxmSyy = Sxx - Syy

    a11 = SxxpSyy + Szz - mxEigenV
    a12 = SyzmSzy
    a13 = -SxzmSzx
    a14 = SxymSyx
    a21 = SyzmSzy
    a22 = SxxmSyy - Szz - mxEigenV
    a23 = SxypSyx
    a24 = SxzpSzx
    a31 = a13
    a32 = a23
    a33 = Syy - Sxx - Szz - mxEigenV
    a34 = SyzpSzy
    a41 = a14
    a42 = a24
    a43 = a34
    a44 = Szz - SxxpSyy - mxEigenV

    a3344_4334 = a33 * a44 - a43 * a34
    a3244_4234 = a32 * a44 - a42 * a34
    a3243_4233 = a32 * a43 - a42 * a33
    a3143_4133 = a31 * a43 - a41 * a33
    a3144_4134 = a31 * a44 - a41 * a34
    a3142_4132 = a31 * a42 - a41 * a32

    q1 =  a22*a3344_4334 - a23*a3244_4234 + a24*a3243_4233
    q2 = -a21*a3344_4334 + a23*a3144_4134 - a24*a3143_4133
    q3 =  a21*a3244_4234 - a22*a3144_4134 + a24*a3142_4132
    q4 = -a21*a3243_4233 + a22*a3143_4133 - a23*a3142_4132

    qsqr = q1*q1 + q2*q2 + q3*q3 + q4*q4

    if qsqr < _EVECPREC:
        q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133
        q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132
        qsqr = q1*q1 + q2*q2 + q3*q3 + q4*q4

        if qsqr < _EVECPREC:
            a1324_1423 = a13*a24 - a14*a23
            a1224_1422 = a12*a24 - a14*a22
            a1223_1322 = a12*a23 - a13*a22
            a1124_1421 = a11*a24 - a14*a21
            a1123_1321 = a11*a23 - a13*a21
            a1122_1221 = a11*a22 - a12*a21

            q1 =  a42*a1324_1423 - a43*a1224_1422 + a44*a1223_1322
            q2 = -a41*a1324_1423 + a43*a1124_1421 - a44*a1123_1321
            q3 =  a41*a1224_1422 - a42*a1124_1421 + a44*a1122_1221
            q4 = -a41*a1223_1322 + a42*a1123_1321 - a43*a1122_1221
            qsqr = q1*q1 + q2*q2 + q3*q3 + q4*q4

            if qsqr < _EVECPREC:
                q1 =  a32*a1324_1423 - a33*a1224_1422 + a34*a1223_1322
                q2 = -a31*a1324_1423 + a33*a1124_1421 - a34*a1123_1321
                q3 =  a31*a1224_1422 - a32*a1124_1421 + a34*a1122_1221
                q4 = -a31*a1223_1322 + a32*a1123_1321 - a33*a1122_1221
                qsqr = q1*q1 + q2*q2 + q3*q3 + q4*q4

                if qsqr < _EVECPREC:
                    rot = np.array([1., 0., 0., 0., 1., 0., 0., 0., 1.])
                    return rmsd, rot

    normq = np.sqrt(qsqr)
    q1 /= normq; q2 /= normq; q3 /= normq; q4 /= normq

    a2 = q1*q1; x2 = q2*q2; y2 = q3*q3; z2 = q4*q4
    xy = q2*q3; az = q1*q4; zx = q4*q2; ay = q1*q3; yz = q3*q4; ax = q1*q2

    rot = np.array([
        a2 + x2 - y2 - z2,  2*(xy + az),  2*(zx - ay),
        2*(xy - az),  a2 - x2 + y2 - z2,  2*(yz + ax),
        2*(zx + ay),  2*(yz - ax),  a2 - x2 - y2 + z2
    ])
    return rmsd, rot


def calc_rmsd(coords1, coords2):
    """
    Calculate the minimum RMSD between two sets of coordinates.

    Parameters
    ----------
    coords1 : array_like, shape (N, 3)
        Reference coordinates.
    coords2 : array_like, shape (N, 3)
        Target coordinates.

    Returns
    -------
    rmsd : float
    """
    c1 = np.asarray(coords1, dtype=np.float64)
    c2 = np.asarray(coords2, dtype=np.float64)
    c1, _ = _center_coords(c1)
    c2, _ = _center_coords(c2)
    A, E0 = _inner_product(c1, c2)
    return _fast_calc_rmsd(A, len(c1), E0)


def calc_rotation(coords1, coords2):
    """
    Calculate RMSD and 4x4 homogeneous rotation-translation matrix.

    The returned matrix transforms coords2 onto coords1:
        coords2_transformed = coords2 @ mat[:3,:3].T + mat[:3,3]

    Parameters
    ----------
    coords1 : array_like, shape (N, 3)
        Reference coordinates.
    coords2 : array_like, shape (N, 3)
        Target coordinates.

    Returns
    -------
    rmsd : float
    mat : ndarray, shape (4, 4)
        Homogeneous rotation-translation matrix.
    """
    c1 = np.asarray(coords1, dtype=np.float64)
    c2 = np.asarray(coords2, dtype=np.float64)
    center1 = c1.mean(axis=0)
    center2 = c2.mean(axis=0)
    c1c = c1 - center1
    c2c = c2 - center2
    A, E0 = _inner_product(c1c, c2c)
    rmsd, rot = _fast_calc_rmsd_and_rotation(A, len(c1), E0)
    rotmat = rot.reshape(3, 3)

    # Build 4x4: mat = TR1 @ ROT @ TR2
    # TR1 translates to center1, ROT rotates, TR2 translates from center2 to origin
    mat = np.eye(4)
    mat[:3, :3] = rotmat
    mat[:3, 3] = center1 - rotmat @ center2
    return rmsd, mat


def superimpose(coords1, coords2):
    """
    Superimpose coords2 onto coords1 and return the transformed coordinates.

    Parameters
    ----------
    coords1 : array_like, shape (N, 3)
        Reference coordinates.
    coords2 : array_like, shape (N, 3)
        Target coordinates to be transformed.

    Returns
    -------
    rmsd : float
    coords2_transformed : ndarray, shape (N, 3)
        coords2 after optimal rotation and translation onto coords1.
    """
    c1 = np.asarray(coords1, dtype=np.float64)
    c2 = np.asarray(coords2, dtype=np.float64)
    center1 = c1.mean(axis=0)
    center2 = c2.mean(axis=0)
    c1c = c1 - center1
    c2c = c2 - center2
    A, E0 = _inner_product(c1c, c2c)
    rmsd, rot = _fast_calc_rmsd_and_rotation(A, len(c1), E0)
    rotmat = rot.reshape(3, 3)

    coords2_tr = (c2 - center2) @ rotmat.T + center1
    return rmsd, coords2_tr
