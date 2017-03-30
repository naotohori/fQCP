#define PREC 8

subroutine InnerProduct(nlen, coords1, coords2, A, E0, weight)
    
   implicit none 

   integer,    intent(in) :: nlen
   real(PREC), intent(in) :: coords1(3,nlen)
   real(PREC), intent(in) :: coords2(3,nlen)
   real(PREC), intent(out) :: A(9)
   real(PREC), intent(out) :: E0
   real(PREC), intent(in), optional :: weight(nlen)

   integer :: i
   real(PREC) :: G

   A(:) = 0.0e0
   G = 0.0e0

   if (present(weight)) then
      do i = 1, nlen
         G = G + dot_product(coords1(1:3,i), coords1(1:3,i)) * weight(i) &
               + dot_product(coords2(1:3,i), coords2(1:3,i)) * weight(i)
            
         A(1:3) = A(1:3) + coords1(1,i) * coords2(1:3,i)
         A(4:6) = A(4:6) + coords1(2,i) * coords2(1:3,i)
         A(7:9) = A(7:9) + coords1(3,i) * coords2(1:3,i)
      enddo
   else
      do i = 1, nlen
         G = G + dot_product( coords1(1:3,i), coords1(1:3,i) ) &
               + dot_product( coords2(1:3,i), coords2(1:3,i) )

         A(1:3) = A(1:3) + coords1(1,i) * coords2(1:3,i)
         A(4:6) = A(4:6) + coords1(2,i) * coords2(1:3,i)
         A(7:9) = A(7:9) + coords1(3,i) * coords2(1:3,i)
      enddo
   endif
   
   E0 = G * 0.5
endsubroutine InnerProduct

subroutine CenterCoords(nlen, coords, weight)
    
   implicit none 

   integer, intent(in) :: nlen
   real(PREC), intent(inout) :: coords(3,nlen)
   real(PREC), intent(in), optional :: weight(nlen)
   
   integer :: i
   real(PREC) :: s(3), wsum

   s(:) = 0.0

   if (present(weight)) then
      wsum = 0.0
      do i = 1, nlen
         s(:) = s(:) + weight(i) * coords(:,i)
         wsum = wsum + weight(i)
      enddo

      s(:) = s(:) / wsum
   else
      do i = 1, nlen
         s(:) = s(:) + coords(:,i)
      enddo

      s(:) = s(:) / real(nlen, kind=PREC)
   endif

   do i = 1, nlen
      coords(:,i) = coords(:,i) - s(:)
   enddo
endsubroutine CenterCoords

subroutine InnerProductNoWeight(nlen, coords1, coords2, A, E0)
    
   implicit none 

   integer,    intent(in) :: nlen
   real(PREC), intent(in) :: coords1(3,nlen)
   real(PREC), intent(in) :: coords2(3,nlen)
   real(PREC), intent(out) :: A(9)
   real(PREC), intent(out) :: E0

   integer :: i
   real(PREC) :: G

   A(:) = 0.0e0
   G = 0.0e0

   do i = 1, nlen
      G = G + dot_product( coords1(1:3,i), coords1(1:3,i) ) &
            + dot_product( coords2(1:3,i), coords2(1:3,i) )

      A(1:3) = A(1:3) + coords1(1,i) * coords2(1:3,i)
      A(4:6) = A(4:6) + coords1(2,i) * coords2(1:3,i)
      A(7:9) = A(7:9) + coords1(3,i) * coords2(1:3,i)
   enddo
   
   E0 = G * 0.5
endsubroutine InnerProductNoWeight

subroutine CenterCoordsNoWeight(nlen, coords)
    
   implicit none 

   integer, intent(in) :: nlen
   real(PREC), intent(inout) :: coords(3,nlen)
   
   integer :: i
   real(PREC) :: s(3)

   s(:) = 0.0

   do i = 1, nlen
      s(:) = s(:) + coords(:,i)
   enddo

   s(:) = s(:) / real(nlen, kind=PREC)

   do i = 1, nlen
      coords(:,i) = coords(:,i) - s(:)
   enddo
endsubroutine CenterCoordsNoWeight
