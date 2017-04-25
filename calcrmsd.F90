real(8) function CalcRMSD(nlen, coords1, coords2)
    
   implicit none 

   integer, intent(in)  :: nlen
   real(8), intent(inout) :: coords1(3,nlen)
   real(8), intent(inout) :: coords2(3,nlen)
   real(8) FastCalcRMSD

   real(8) :: A(9), E0

   call CenterCoordsNoWeight(nlen, coords1)
   call CenterCoordsNoWeight(nlen, coords2)

   call InnerProductNoWeight(nlen, coords1, coords2, A, E0)

   ! calculate the RMSD & rotational matrix
   CalcRMSD =  FastCalcRMSD(A, nlen, E0)
   return

endfunction CalcRMSD
