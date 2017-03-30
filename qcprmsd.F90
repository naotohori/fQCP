! Fotran version of Quaternion Characteristic Polynomial (QCP) 
! Copyright (c) 2016 Naoto Hori
!
! The original code written in c is distributed at http://theobald.brandeis.edu/qcp
! See further information in README.md
!
! Following is the original notice by the original authors.
!
!/*******************************************************************************
! *  -/_|:|_|_\- 
! *
! *  File:           qcprot.c
! *  Version:        1.4
! *
! *  Function:       Rapid calculation of the least-squares rotation using a 
! *                  quaternion-based characteristic polynomial and 
! *                  a cofactor matrix
! *
! *  Author(s):      Douglas L. Theobald
! *                  Department of Biochemistry
! *                  MS 009
! *                  Brandeis University
! *                  415 South St
! *                  Waltham, MA  02453
! *                  USA
! *
! *                  dtheobald@brandeis.edu
! *                  
! *                  Pu Liu
! *                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
! *                  665 Stockton Drive
! *                  Exton, PA  19341
! *                  USA
! *
! *                  pliu24@its.jnj.com
! * 
! *
! *    If you use this QCP rotation calculation method in a publication, please
! *    reference:
! *
! *      Douglas L. Theobald (2005)
! *      "Rapid calculation of RMSD using a quaternion-based characteristic
! *      polynomial."
! *      Acta Crystallographica A 61(4):478-480.
! *
! *      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
! *      "Fast determination of the optimal rotational matrix for macromolecular 
! *      superpositions."
! *      Journal of Computational Chemistry 31(7):1561-1563.
! *
! *
! *  Copyright (c) 2009-2013 Pu Liu and Douglas L. Theobald
! *  All rights reserved.
! *
! *  Redistribution and use in source and binary forms, with or without modification, are permitted
! *  provided that the following conditions are met:
! *
! *  * Redistributions of source code must retain the above copyright notice, this list of
! *    conditions and the following disclaimer.
! *  * Redistributions in binary form must reproduce the above copyright notice, this list
! *    of conditions and the following disclaimer in the documentation and/or other materials
! *    provided with the distribution.
! *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to
! *    endorse or promote products derived from this software without specific prior written
! *    permission.
! *
! *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
! *
! *  Source:         started anew.
! *
! *  Change History:
! *    2009/04/13      Started source
! *    2010/03/28      Modified FastCalcRMSDAndRotation() to handle tiny qsqr
! *                    If trying all rows of the adjoint still gives too small
! *                    qsqr, then just return identity matrix. (DLT)
! *    2010/06/30      Fixed prob in assigning A[9] = 0 in InnerProduct()
! *                    invalid mem access
! *    2011/02/21      Made CenterCoords use weights
! *    2011/05/02      Finally changed CenterCoords declaration in qcprot.h
! *                    Also changed some functions to static
! *    2011/07/08      put in fabs() to fix taking sqrt of small neg numbers, fp error
! *    2012/07/26      minor changes to comments and main.c, more info (v.1.4)
! *  
! ******************************************************************************/

#define PREC 8

real(PREC) function FastCalcRMSD(A, nlen, E0)
    
   implicit none 

   real(PREC), intent(in)  :: A(9)
   integer,    intent(in)  :: nlen
   real(PREC), intent(in)  :: E0

   real(PREC) :: Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz
   real(PREC) :: Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2, &
                 SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2, &
                 SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy, &
                 SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy
   real(PREC) :: C(3)
   integer :: i
   real(PREC) :: mxEigenV
   real(PREC) :: oldg
   real(PREC) :: b, aa, delta
   real(PREC) :: x2
   real(PREC), parameter :: evalprec = 1.0e-11

   Sxx = A(1)
   Sxy = A(2)
   Sxz = A(3)
   Syx = A(4)
   Syy = A(5)
   Syz = A(6)
   Szx = A(7)
   Szy = A(8)
   Szz = A(9)

   Sxx2 = Sxx * Sxx
   Syy2 = Syy * Syy
   Szz2 = Szz * Szz

   Sxy2 = Sxy * Sxy
   Syz2 = Syz * Syz
   Sxz2 = Sxz * Sxz

   Syx2 = Syx * Syx
   Szy2 = Szy * Szy
   Szx2 = Szx * Szx

   SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz)
   Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

   C(3) = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
   C(2) =  8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz)

   SxzpSzx = Sxz + Szx
   SyzpSzy = Syz + Szy
   SxypSyx = Sxy + Syx
   SyzmSzy = Syz - Szy
   SxzmSzx = Sxz - Szx
   SxymSyx = Sxy - Syx
   SxxpSyy = Sxx + Syy
   SxxmSyy = Sxx - Syy
   Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

   C(1) = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 &
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) &
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz)) &
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz)) &
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz)) &
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz))

   ! Newton-Raphson
   mxEigenV = E0
   do i = 1, 50
      oldg = mxEigenV
      x2 = mxEigenV*mxEigenV
      b = (x2 + C(3))*mxEigenV
      aa = b + C(2)
      delta = ((aa*mxEigenV + C(1))/(2.0*x2*mxEigenV + b + aa))
      mxEigenV = mxEigenV - delta
      ! write(*,*) "diff[",i,"]:",mxEigenV-oldg, evalprec*mxEigenV, mxEigenV
      if (abs(mxEigenV - oldg) < abs(evalprec*mxEigenV)) then
         exit
      endif
   enddo

   if (i == 50) then
      !write(*,*) "More than",i,"iterations needed!"
      FastCalcRMSD = -1.0
      return
   endif

   ! the abs() is to guard against extremely small, but *negative* numbers due to floating point error
   FastCalcRMSD = sqrt(abs(2.0 * (E0 - mxEigenV)/real(nlen,kind=PREC)))
   return
endfunction FastCalcRMSD

