subroutine Superimpose(nlen, coords1, coords2, rmsd, rot, weight)
    
#define PREC 8

   implicit none 

   interface
      subroutine CenterCoords(nlen, coords, weight)
         integer, intent(in) :: nlen
         real(PREC), intent(inout) :: coords(3,nlen)
         real(PREC), intent(in), optional :: weight(nlen)
      endsubroutine CenterCoords
      subroutine InnerProduct(nlen, coords1, coords2, A, E0, weight)
         integer,    intent(in) :: nlen
         real(PREC), intent(in) :: coords1(3,nlen)
         real(PREC), intent(in) :: coords2(3,nlen)
         real(PREC), intent(out) :: A(9)
         real(PREC), intent(out) :: E0
         real(PREC), intent(in), optional :: weight(nlen)
      endsubroutine InnerProduct
      subroutine FastCalcRMSDAndRotation(A, minScore, nlen, E0, rmsd, rot, ireturn)
         real(PREC), intent(in)  :: A(9)
         integer,    intent(in)  :: minScore
         integer,    intent(in)  :: nlen
         real(PREC), intent(in)  :: E0
         real(PREC), intent(out) :: rmsd
         real(PREC), intent(out) :: rot(9)
         integer,    intent(out) :: ireturn
      endsubroutine FastCalcRMSDAndRotation
   endinterface

   integer, intent(in)  :: nlen
   real(8), intent(inout) :: coords1(3,nlen)
   real(8), intent(inout) :: coords2(3,nlen)
   real(8), intent(out) :: rmsd
   real(8), intent(out) :: rot(9)
   real(8), intent(in), optional :: weight(nlen)

   integer :: ireturn
   real(8) :: A(9), E0

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
   call FastCalcRMSDAndRotation(A, -1, nlen, E0, rmsd, rot, ireturn)

endsubroutine Superimpose
