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
    real(c_double), pointer :: X(:,:), KS(:), Ntrl_r, Nthread_r
    integer(c_int16_t), pointer :: Y(:)
    mwPointer :: mxKS
    mwSize, parameter :: One=1
    mwSize :: Nx, Ntrl, Nthread

    if( nrhs /= 4 ) then
        call mexErrMsgTxt("This function takes 5 inputs")
    endif
    if( nlhs > 1) then
        call mexErrMsgTxt("This function returns 1 output")
    endif

    call fpGetPr( X,  prhs(1) )
    call fpGetPr( Y,  prhs(2) )
    call fpGetPr( Ntrl_r,  prhs(3) )
    call fpGetPr( Nthread_r,  prhs(4) )

    if( (.not.associated(X)) .or. &
        (.not.associated(Y)) .or. (.not.associated(Ntrl_r)) .or. &
        (.not.associated(Nthread_r)) ) then
        call mexErrMsgTxt("Problem with inputs: check types and dimensions")
    endif

    Ntrl = Ntrl_r
    Nthread = Nthread_r
    if( (size(X,1)/=Ntrl) .or. (size(Y)/=Ntrl) ) then
        call mexErrMsgTxt("Number of trials does not match data")
    endif
    Nx = size(X,2)

    ! allocate Matlab memory
    mxKS = mxCreateNumericMatrix(One, Nx, mxDOUBLE_CLASS, mxREAL)
    call fpGetPr( KS, mxKS )
    plhs(1) = mxKS

    call calc_kstest_slice_omp(X, Nx, Y, Ntrl, Nthread, KS)


end subroutine mexFunction


subroutine calc_kstest_slice_omp(X, Nx, Y, Ntrl, Nthread, KS)
    use iso_c_binding
    use MatlabAPImx
    use MatlabAPImex
    use lib_array
    implicit none

    ! ARG
    mwSize, intent(in) :: Nx, Ntrl, Nthread
    real(c_double), intent(in) :: X(Ntrl,Nx)
    integer(c_int16_t), intent(in) :: Y(Ntrl)
    real(c_double), intent(out) :: KS(Nx)
    ! LOC
    mwSize, parameter :: One=1
    integer :: idx(Ntrl), xi, ti, N0, N1, i
    real(c_double) :: cX(Ntrl), cdf0, cdf1, N0inv, N1inv, ksdiff(Ntrl)
    character(len=255) :: dbgmsg

    N0 = 0
    N1 = 0
    do ti=1,Ntrl
      if( Y(ti)==0 ) then
        N0 = N0 + 1
      else
        N1 = N1 + 1
      endif
    end do
    N0inv = 1.0_8 / real(N0,8)
    N1inv = 1.0_8 / real(N1,8)

 
    ! loop over X variables
    ! NB: OMP doesn't work with allocatable arrays
    ! NB: default(private) shared(...list...) didn't work
    !$omp parallel do &
    !$omp num_threads(Nthread) &
    !$omp default(shared) &
    !$omp private(xi,idx,i,cX,ti,cdf0,cdf1,ksdiff)
    do xi=1,Nx
      ! argsort data
      idx = (/ (i, i=1,Ntrl) /)
      cX = X(:,xi)
      call quicksort_index(cX, idx)
      cdf0 = 0.0_8
      cdf1 = 0.0_8

      if( Y(idx(1))==0 ) then
        cdf0 = N0inv
      else
        cdf1 = N1inv
      endif
      ksdiff(1) = abs(cdf0 - cdf1)

      do ti=2,Ntrl
        if( Y(idx(ti))==0 ) then
          cdf0 = cdf0 + N0inv
          cdf1 = cdf1
        else
          cdf0 = cdf0
          cdf1 = cdf1 + N1inv
        endif
        ksdiff(ti) = abs(cdf0 - cdf1)
      end do

      KS(xi) = maxval(ksdiff)
    end do
    !$omp end parallel do

end subroutine
