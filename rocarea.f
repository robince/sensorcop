#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use iso_c_binding
    use MatlabAPImx
    use MatlabAPImex
    implicit none
! ARG
    mwPointer :: plhs(*), prhs(*)
    integer :: nlhs, nrhs
! LOC
    real(c_double), pointer :: p(:), Az
    integer(c_int16_t), pointer :: label(:)
    mwPointer :: mxAz
    mwSize, parameter :: One=1
    mwSize :: Ntrl

    if( nrhs /= 2 ) then
        call mexErrMsgTxt("This function takes 2 inputs")
    endif
    if( nlhs > 1) then
        call mexErrMsgTxt("This function returns 1 output")
    endif

    call fpGetPr( p,  prhs(1) )
    call fpGetPr( label,  prhs(2) )

    if( (.not.associated(p)) .or. (.not.associated(label)) ) then
        call mexErrMsgTxt("Problem with inputs: check types and dimensions")
    endif

    Ntrl = size(p)
    if( (size(label)/=Ntrl) ) then
        call mexErrMsgTxt("Number of trials does not match")
    endif

    ! allocate Matlab memory
    mxAz = mxCreateNumericMatrix(One, One, mxDOUBLE_CLASS, mxREAL)
    call fpGetPr( Az, mxAz )
    plhs(1) = mxAz

    call dorocarea(p, label, Ntrl, Az)

end subroutine mexFunction

subroutine dorocarea(p, label, Ntrl, Az)
    use iso_c_binding
    use lib_array
    implicit none
    ! ARG
    mwSize, intent(in) :: Ntrl
    real(c_double), intent(in) :: p(Ntrl)
    integer(c_int16_t), intent(in) :: label(Ntrl)
    real(c_double), intent(out) :: Az
    ! LOC
    integer :: index(Ntrl), i, Np, Nn
    real(c_double) :: cp(Ntrl), tp, denNp, denNn

    Nn = 0
    Np = 0
    do i=1,Ntrl
      if( label(i) == 0 ) then 
        Nn = Nn+1
      else 
        Np = Np+1
      endif
    end do

    Az = 0.0

    ! argsort data
    index = (/ (i, i=1,Ntrl) /)
    cp = -p
    call quicksort_index(cp,index)

    Az = 0.0
    tp = 0.0
    denNp = 1.0_8 / real(Np,8)
    denNn = 1.0_8 / real(Nn,8)
    do i=1,Ntrl
      if( label(index(i))==1 ) then
        tp = tp + label(index(i))*denNp
      else
        ! label(index(i)) == 0
        Az = Az + tp*denNn
      endif
    end do
end subroutine dorocarea

