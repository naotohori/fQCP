subroutine Superimpose(nlen, coords1_io, coords2_io, rmsd)
    
   implicit none 

   integer, intent(in)  :: nlen
   real(8), intent(inout) :: coords1_io(3,nlen)
   real(8), intent(inout) :: coords2_io(3,nlen)
   real(8), intent(out) :: rmsd

   integer :: i
   integer :: ireturn
   real(8) :: A(9), E0
   real(8) :: rot(9), rotmat(3,3)
   real(8) :: coords1(3,nlen), coords2(3,nlen)
   real(8) :: center1(3), center2(3)

   coords1(:,:) = coords1_io(:,:)
   coords2(:,:) = coords2_io(:,:)

   call CenterCoordsNoWeight(nlen, coords1)
   call CenterCoordsNoWeight(nlen, coords2)

   call InnerProductNoWeight(nlen, coords1, coords2, A, E0)

   ! calculate the RMSD & rotational matrix
   call FastCalcRMSDAndRotation(A, -1, nlen, E0, rmsd, rot, ireturn)


   ! Apply translation to coords2
   center1(:) = 0.0
   center2(:) = 0.0
   do i = 1, nlen
      center1(:) = center1(:) + coords1_io(:,i)
      center2(:) = center2(:) + coords2_io(:,i)
   enddo
   center1(:) = center1(:) / real(nlen)
   center2(:) = center2(:) / real(nlen)

   ! Apply rotation to coords2
   !rotmat(1,1) = rot(1)
   !rotmat(1,2) = rot(2)
   !rotmat(1,3) = rot(3)
   !rotmat(2,1) = rot(4)
   !rotmat(2,2) = rot(5)
   !rotmat(2,3) = rot(6)
   !rotmat(3,1) = rot(7)
   !rotmat(3,2) = rot(8)
   !rotmat(3,3) = rot(9)
   rotmat(1,1) = rot(1)
   rotmat(2,1) = rot(2)
   rotmat(3,1) = rot(3)
   rotmat(1,2) = rot(4)
   rotmat(2,2) = rot(5)
   rotmat(3,2) = rot(6)
   rotmat(1,3) = rot(7)
   rotmat(2,3) = rot(8)
   rotmat(3,3) = rot(9)

   do i = 1, nlen
      coords2_io(:,i) = coords2_io(:,i) - center2(:)
      coords2_io(:,i) = matmul(rotmat, coords2_io(:,i))
      coords2_io(:,i) = coords2_io(:,i) + center1(:)
   enddo

endsubroutine Superimpose
