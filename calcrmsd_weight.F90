real(8) function CalcRMSD(nlen, coords1, coords2, weight)
    
   implicit none 

   integer, intent(in)  :: nlen
   real(8), intent(inout) :: coords1(3,nlen)
   real(8), intent(inout) :: coords2(3,nlen)
   real(8), intent(in), optional :: weight(nlen)

   real(8) :: A(9), E0

   interface
      subroutine CenterCoords(nlen, coords, weight)
         integer, intent(in) :: nlen
         real(8), intent(inout) :: coords(3,nlen)
         real(8), intent(in), optional :: weight(nlen)
      endsubroutine CenterCoords
      subroutine InnerProduct(nlen, coords1, coords2, A, E0, weight)
         integer,    intent(in) :: nlen
         real(8), intent(in) :: coords1(3,nlen)
         real(8), intent(in) :: coords2(3,nlen)
         real(8), intent(out) :: A(9)
         real(8), intent(out) :: E0
         real(8), intent(in), optional :: weight(nlen)
      endsubroutine InnerProduct
      real(8) function FastCalcRMSD(A, nlen, E0)
         real(8), intent(in)  :: A(9)
         integer,    intent(in)  :: nlen
         real(8), intent(in)  :: E0
      endfunction FastCalcRMSD
   endinterface

   if (present(weight)) then
      ! center the structures -- if precentered you can omit this step
      call CenterCoords(nlen, coords1, weight)
      call CenterCoords(nlen, coords2, weight)
       
      ! calculate the (weighted) inner product of two structures
      call InnerProduct(nlen, coords1, coords2, A, E0, weight)
   else
      call CenterCoords(nlen, coords1)
      call CenterCoords(nlen, coords2)

      call InnerProduct(nlen, coords1, coords2, A, E0)
   endif

   ! calculate the RMSD & rotational matrix
   CalcRMSD =  FastCalcRMSD(A, nlen, E0)
   return

endfunction CalcRMSD
