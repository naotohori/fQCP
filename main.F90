! /* Sample code to use the routine for fast RMSD & rotational matrix calculation.
!    Note that we superposition frag_b onto frag_a.  
!    For the example provided below, the minimum least-squares RMSD for the two
!    7-atom fragments should be 0.719106 A.
!
!    The rotation quaterion should be:
!
!    -8.620063e-01   3.435505e-01   1.242953e-01  -3.513814e-01 
!
!    And the corresponding 3x3 rotation matrix is:
!
!    [     0.72216358     0.69118937    -0.02714790 ]
!    [    -0.52038257     0.51700833    -0.67963547 ]
!    [    -0.45572112     0.50493528     0.73304748 ]
!*/

#define PREC 8

program testQCP

   implicit none

   interface
      subroutine CalcRMSDRotationalMatrix(nlen, coords1, coords2, rmsd, rot, weight)
         integer, intent(in)  :: nlen
         real(PREC), intent(inout) :: coords1(3,nlen)
         real(PREC), intent(inout) :: coords2(3,nlen)
         real(PREC), intent(out) :: rmsd
         real(PREC), intent(out) :: rot(9)
         real(PREC), intent(in), optional :: weight(nlen)
      endsubroutine CalcRMSDRotationalMatrix
   endinterface
   integer, parameter :: nlen = 7

   integer :: i
   real(PREC) :: rmsd, x, y, z, euc_dist
   real(PREC) :: frag_a(3,nlen), frag_b(3,nlen)
   real(PREC) :: rotmat(9), d(3)

   frag_a(1,1) =  -2.803
   frag_a(2,1) = -15.373
   frag_a(3,1) =  24.556
   frag_a(1,2) =   0.893
   frag_a(2,2) = -16.062
   frag_a(3,2) =  25.147
   frag_a(1,3) =   1.368
   frag_a(2,3) = -12.371
   frag_a(3,3) =  25.885
   frag_a(1,4) =  -1.651
   frag_a(2,4) = -12.153
   frag_a(3,4) =  28.177
   frag_a(1,5) =  -0.440
   frag_a(2,5) = -15.218
   frag_a(3,5) =  30.068
   frag_a(1,6) =   2.551
   frag_a(2,6) = -13.273
   frag_a(3,6) =  31.372
   frag_a(1,7) =   0.105
   frag_a(2,7) = -11.330
   frag_a(3,7) =  33.567

   frag_b(1,1) = -14.739
   frag_b(2,1) = -18.673
   frag_b(3,1) =  15.040
   frag_b(1,2) = -12.473
   frag_b(2,2) = -15.810
   frag_b(3,2) =  16.074
   frag_b(1,3) = -14.802
   frag_b(2,3) = -13.307
   frag_b(3,3) =  14.408
   frag_b(1,4) = -17.782
   frag_b(2,4) = -14.852
   frag_b(3,4) =  16.171
   frag_b(1,5) = -16.124
   frag_b(2,5) = -14.617
   frag_b(3,5) =  19.584
   frag_b(1,6) = -15.029
   frag_b(2,6) = -11.037
   frag_b(3,6) =  18.902
   frag_b(1,7) = -18.577
   frag_b(2,7) = -10.001
   frag_b(3,7) =  17.996

   write(*,*) ''
   write(*,*) "Coords before centering:"
   write(*,*) ''
   do i = 1, nlen
      write(*,'(f8.3 f8.3 f8.3)') frag_a(1,i), frag_a(2,i), frag_a(3,i)
   enddo
   write(*,*) ''
   do i = 1, nlen
      write(*,'(f8.3 f8.3 f8.3)') frag_b(1,i), frag_b(2,i), frag_b(3,i)
   enddo
   write(*,*) ''

   call CalcRMSDRotationalMatrix(nlen, frag_a, frag_b, rmsd, rotmat)

   write(*,*) ''
   write(*,*) "Coords after centering:"
   write(*,*) ''
   do i = 1, nlen
      write(*,'(f8.3 f8.3 f8.3)') frag_a(1,i), frag_a(2,i), frag_a(3,i)
   enddo
   write(*,*) ''
   do i = 1, nlen
      write(*,'(f8.3 f8.3 f8.3)') frag_b(1,i), frag_b(2,i), frag_b(3,i)
   enddo
   write(*,*) ''

   write(*,*) ''
   write(*,*) "QCP rmsd:", rmsd
   write(*,*) ''

   write(*,*) ''
   write(*,*) "QCP Rotation matrix:"
   write(*,*) ''

   write(*,'(f14.8 f14.8 f14.8)') rotmat(1), rotmat(2), rotmat(3)
   write(*,'(f14.8 f14.8 f14.8)') rotmat(4), rotmat(5), rotmat(6)
   write(*,'(f14.8 f14.8 f14.8)') rotmat(7), rotmat(8), rotmat(9)
   write(*,*) ''

   ! apply rotation matrix
   do i = 1, nlen
      x = dot_product(rotmat(1:3), frag_b(1:3,i))
      y = dot_product(rotmat(4:6), frag_b(1:3,i))
      z = dot_product(rotmat(7:9), frag_b(1:3,i))
      
      frag_b(1,i) = x
      frag_b(2,i) = y
      frag_b(3,i) = z
   enddo

   ! calculate euclidean distance
   euc_dist = 0.0

   do i = 1, nlen
      d(1:3) = frag_a(1:3,i) - frag_b(1:3,i)
      euc_dist = euc_dist + dot_product(d,d)
   enddo

   write(*,*) ''
   write(*,*) "Coords 2 after rotation:"
   write(*,*) ''
   do i = 1, nlen
      write(*,'(f8.3 f8.3 f8.3)') frag_b(1,i), frag_b(2,i), frag_b(3,i)
   enddo
   write(*,*) ''

   write(*,*) ''
   write(*,*) "Explicit RMSD calculated from transformed coords: ", sqrt(euc_dist / real(nlen,kind=PREC))
   write(*,*) ''
   write(*,*) ''

endprogram testQCP
