!subroutine CalcRotation(nlen, coords1_in, coords2_in, rmsd, rotmat, weight)
subroutine CalcRotation(nlen, coords1_in, coords2_in, rmsd, mat)

   implicit none 

!   interface
!      subroutine CenterCoords(nlen, coords, weight)
!         integer, intent(in) :: nlen
!         real(8), intent(inout) :: coords(3,nlen)
!         real(8), intent(in), optional :: weight(nlen)
!      endsubroutine CenterCoords
!      subroutine InnerProduct(nlen, coords1, coords2, A, E0, weight)
!         integer,    intent(in) :: nlen
!         real(8), intent(in) :: coords1(3,nlen)
!         real(8), intent(in) :: coords2(3,nlen)
!         real(8), intent(out) :: A(9)
!         real(8), intent(out) :: E0
!         real(8), intent(in), optional :: weight(nlen)
!      endsubroutine InnerProduct
!      subroutine FastCalcRMSDAndRotation(A, minScore, nlen, E0, rmsd, rot, ireturn)
!         real(8), intent(in)  :: A(9)
!         integer,    intent(in)  :: minScore
!         integer,    intent(in)  :: nlen
!         real(8), intent(in)  :: E0
!         real(8), intent(out) :: rmsd
!         real(8), intent(out) :: rot(9)
!         integer,    intent(out) :: ireturn
!      endsubroutine FastCalcRMSDAndRotation
!   endinterface

   integer, intent(in)  :: nlen
   real(8), intent(in) :: coords1_in(3,nlen)
   real(8), intent(in) :: coords2_in(3,nlen)
   real(8), intent(out) :: rmsd
   real(8), intent(out) :: mat(4,4)
!   real(8), intent(in), optional :: weight(nlen)

   integer :: ireturn, i
   real(8) :: rot(9)
   real(8) :: A(9), E0
   real(8) :: coords1(3,nlen)
   real(8) :: coords2(3,nlen)
   real(8) :: center1(3), center2(3)
   real(8) :: rotmat(4,4)

   coords1(:,:) = coords1_in(:,:)
   coords2(:,:) = coords2_in(:,:)

!   if (present(weight)) then
!      ! center the structures -- if precentered you can omit this step
!      call CenterCoords(nlen, coords1, weight)
!      call CenterCoords(nlen, coords2, weight)
!       
!      ! calculate the (weighted) inner product of two structures
!      call InnerProduct(nlen, coords1, coords2, A, E0, weight)
!   else
      call CenterCoordsNoWeight(nlen, coords1)
      call CenterCoordsNoWeight(nlen, coords2)

      call InnerProductNoWeight(nlen, coords1, coords2, A, E0)
!   endif

   ! calculate the RMSD & rotational matrix
   call FastCalcRMSDAndRotation(A, -1, nlen, E0, rmsd, rot, ireturn)


   ! For translation
   center1(:) = 0.0
   center2(:) = 0.0
   do i = 1, nlen
      center1(:) = center1(:) + coords1_in(:,i)
      center2(:) = center2(:) + coords2_in(:,i)
   enddo
   center1(:) = center1(:) / real(nlen)
   center2(:) = center2(:) / real(nlen)


   ! Calculate final rotation-translation matrix
   ! [Trans Origin --> R1] [Rotation] [Trans R2 --> Origin] [Target]
   ! [       TR1         ] [  rot   ] [        TR2        ] [Target]
   ! ===> [mat] [Target]

   ! [ TR1 ]
   mat(:,:) = 0.0
   do i = 1, 4
      mat(i,i) = 1.0
   enddo
   !write(*,*) 'mat 0'
   !do i = 1, 4
   !   write(*,*) mat(i,1), mat(i,2), mat(i,3), mat(i,4) 
   !enddo

   mat(1,4) = center1(1)
   mat(2,4) = center1(2)
   mat(3,4) = center1(3)
   !write(*,*) 'mat 1'
   !do i = 1, 4
   !   write(*,*) mat(i,1), mat(i,2), mat(i,3), mat(i,4) 
   !enddo

   ! [ rot ]
   rotmat(:,:) = 0.0
   do i = 1, 4
      rotmat(i,i) = 1.0
   enddo
   rotmat(1,1) = rot(1)
   rotmat(1,2) = rot(2)
   rotmat(1,3) = rot(3)
   rotmat(2,1) = rot(4)
   rotmat(2,2) = rot(5)
   rotmat(2,3) = rot(6)
   rotmat(3,1) = rot(7)
   rotmat(3,2) = rot(8)
   rotmat(3,3) = rot(9)

   ! Calc [ TR1 ] [ rot ]
   mat = matmul(mat, rotmat)
   !write(*,*) 'mat 2'
   !do i = 1, 4
   !   write(*,*) mat(i,1), mat(i,2), mat(i,3), mat(i,4) 
   !enddo

   ! [ TR2 ]
   rotmat(:,:) = 0.0
   do i = 1, 4
      rotmat(i,i) = 1.0
   enddo
   rotmat(1,4) = - center2(1)
   rotmat(2,4) = - center2(2)
   rotmat(3,4) = - center2(3)

   ! Calc ([ TR1 ] [ rot ]) [ TR2 ]
   mat = matmul(mat, rotmat)
   !write(*,*) 'mat 3'
   !do i = 1, 4
   !   write(*,*) mat(i,1), mat(i,2), mat(i,3), mat(i,4) 
   !enddo

endsubroutine CalcRotation
