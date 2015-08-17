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
    real(c_double), pointer :: X(:,:,:), KS(:), Ntrl_r, Nthread_r
    integer(c_int16_t), pointer :: Y(:)
    mwPointer :: mxKS
    mwSize, parameter :: One=1
    mwSize :: Nx, Ntrl, Nthread, xd

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
    if( (size(X,2)/=2) ) then
      call mexErrMsgTxt("Second dimension of input X must have size 2")
    endif
    Nx = size(X,3)
    xd = 2

    ! allocate Matlab memory
    mxKS = mxCreateNumericMatrix(One, Nx, mxDOUBLE_CLASS, mxREAL)
    call fpGetPr( KS, mxKS )
    plhs(1) = mxKS

    call calc_kstest_slice_omp(X, xd, Nx, Y, Ntrl, Nthread, KS)

end subroutine mexFunction


subroutine calc_kstest_slice_omp(X, xd, Nx, Y, Ntrl, Nthread, KS)
    use iso_c_binding
    use MatlabAPImx
    use MatlabAPImex
    use lib_array
    implicit none

    ! ARG
    mwSize, intent(in) :: Nx, Ntrl, Nthread, xd
    real(c_double), intent(in) :: X(Ntrl,xd,Nx)
    integer(c_int16_t), intent(in) :: Y(Ntrl)
    real(c_double), intent(out) :: KS(Nx)
    ! LOC
    mwSize, parameter :: One=1
    integer :: idx1(Ntrl), idx2(Ntrl), xi, ii, ti, ti1, ti2, N0, N1, i, Y1, Y2
    real(c_double) :: cX1(Ntrl), cX2(Ntrl), cdf0(Ntrl,Ntrl), cdf1(Ntrl,Ntrl), N0inv, N1inv, ksdiff(4)
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
    N0inv = 1.0_8 / (real(N0,8)*real(N0,8))
    N1inv = 1.0_8 / (real(N1,8)*real(N1,8))

 
    ! loop over X variables
    ! NB: OMP doesn't work with allocatable arrays
    ! NB: default(private) shared(...list...) didn't work
    !$omp parallel do &
    !$omp num_threads(Nthread) &
    !$omp default(shared) &
    !$omp private(xi,idx1,idx2,i,cX1,cX2,ti,cdf0,cdf1,ti1,ti2,Y1,Y2,ksdiff)
    do xi=1,Nx
      ! argsort each dimension
      idx1 = (/ (i, i=1,Ntrl) /)
      cX1 = X(:,1,xi)
      call quicksort_index(cX1, idx1)
      idx2 = (/ (i, i=1,Ntrl) /)
      cX2 = X(:,2,xi)
      call quicksort_index(cX2, idx2)

      ! x<X y<Y
      cdf0 = 0.0_8
      cdf1 = 0.0_8
      ! bivariate histogram 
      do ti1=1,Ntrl
        do ti2=1,Ntrl
          Y1 = Y(idx1(ti1))
          Y2 = Y(idx2(ti2))
          if( (Y2==0) .and. (Y1==0) ) then
            cdf0(ti2,ti1) = N0inv
          elseif( (Y2/=0) .and. (Y1/=0) ) then
            cdf1(ti2,ti1) = N1inv
          endif
        end do
      end do
      ! cumsum rows
      do ti1=1,(Ntrl)
        do ti2=1,(Ntrl-1)
          cdf0(ti2+1,ti1) = cdf0(ti2+1,ti1) + cdf0(ti2,ti1)
        end do
      end do
      do ti1=1,(Ntrl)
        do ti2=1,(Ntrl-1)
          cdf1(ti2+1,ti1) = cdf1(ti2+1,ti1) + cdf1(ti2,ti1)
        end do
      end do
      ! cumsum cols
      do ti1=1,(Ntrl-1)
        do ti2=1,(Ntrl)
          cdf0(ti2,ti1+1) = cdf0(ti2,ti1+1) + cdf0(ti2,ti1)
        end do
      end do
      do ti1=1,(Ntrl-1)
        do ti2=1,(Ntrl)
          cdf1(ti2,ti1+1) = cdf1(ti2,ti1+1) + cdf1(ti2,ti1)
        end do
      end do
      ksdiff(1) = maxval(abs(cdf0-cdf1))

      ! x>X y>Y
      cdf0 = 0.0_8
      cdf1 = 0.0_8
      ! bivariate histogram 
      do ti1=1,Ntrl
        do ti2=1,Ntrl
          Y1 = Y(idx1(Ntrl-ti1+1))
          Y2 = Y(idx2(Ntrl-ti2+1))
          if( (Y2==0) .and. (Y1==0) ) then
            cdf0(ti2,ti1) = N0inv
          elseif( (Y2/=0) .and. (Y1/=0) ) then
            cdf1(ti2,ti1) = N1inv
          endif
        end do
      end do
      ! cumsum rows
      do ti1=1,(Ntrl)
        do ti2=1,(Ntrl-1)
          cdf0(ti2+1,ti1) = cdf0(ti2+1,ti1) + cdf0(ti2,ti1)
        end do
      end do
      do ti1=1,(Ntrl)
        do ti2=1,(Ntrl-1)
          cdf1(ti2+1,ti1) = cdf1(ti2+1,ti1) + cdf1(ti2,ti1)
        end do
      end do
      ! cumsum cols
      do ti1=1,(Ntrl-1)
        do ti2=1,(Ntrl)
          cdf0(ti2,ti1+1) = cdf0(ti2,ti1+1) + cdf0(ti2,ti1)
        end do
      end do
      do ti1=1,(Ntrl-1)
        do ti2=1,(Ntrl)
          cdf1(ti2,ti1+1) = cdf1(ti2,ti1+1) + cdf1(ti2,ti1)
        end do
      end do
      ksdiff(2) = maxval(abs(cdf0-cdf1))

      ! x<X y>Y
      cdf0 = 0.0_8
      cdf1 = 0.0_8
      ! bivariate histogram 
      do ti1=1,Ntrl
        do ti2=1,Ntrl
          Y1 = Y(idx1(ti1))
          Y2 = Y(idx2(Ntrl-ti2+1))
          if( (Y2==0) .and. (Y1==0) ) then
            cdf0(ti2,ti1) = N0inv
          elseif( (Y2/=0) .and. (Y1/=0) ) then
            cdf1(ti2,ti1) = N1inv
          endif
        end do
      end do
      ! cumsum rows
      do ti1=1,(Ntrl)
        do ti2=1,(Ntrl-1)
          cdf0(ti2+1,ti1) = cdf0(ti2+1,ti1) + cdf0(ti2,ti1)
        end do
      end do
      do ti1=1,(Ntrl)
        do ti2=1,(Ntrl-1)
          cdf1(ti2+1,ti1) = cdf1(ti2+1,ti1) + cdf1(ti2,ti1)
        end do
      end do
      ! cumsum cols
      do ti1=1,(Ntrl-1)
        do ti2=1,(Ntrl)
          cdf0(ti2,ti1+1) = cdf0(ti2,ti1+1) + cdf0(ti2,ti1)
        end do
      end do
      do ti1=1,(Ntrl-1)
        do ti2=1,(Ntrl)
          cdf1(ti2,ti1+1) = cdf1(ti2,ti1+1) + cdf1(ti2,ti1)
        end do
      end do
      ksdiff(3) = maxval(abs(cdf0-cdf1))

      ! x>X y<Y
      cdf0 = 0.0_8
      cdf1 = 0.0_8
      ! bivariate histogram 
      do ti1=1,Ntrl
        do ti2=1,Ntrl
          Y1 = Y(idx1(Ntrl-ti1+1))
          Y2 = Y(idx2(ti2))
          if( (Y2==0) .and. (Y1==0) ) then
            cdf0(ti2,ti1) = N0inv
          elseif( (Y2/=0) .and. (Y1/=0) ) then
            cdf1(ti2,ti1) = N1inv
          endif
        end do
      end do
      ! cumsum rows
      do ti1=1,(Ntrl)
        do ti2=1,(Ntrl-1)
          cdf0(ti2+1,ti1) = cdf0(ti2+1,ti1) + cdf0(ti2,ti1)
        end do
      end do
      do ti1=1,(Ntrl)
        do ti2=1,(Ntrl-1)
          cdf1(ti2+1,ti1) = cdf1(ti2+1,ti1) + cdf1(ti2,ti1)
        end do
      end do
      ! cumsum cols
      do ti1=1,(Ntrl-1)
        do ti2=1,(Ntrl)
          cdf0(ti2,ti1+1) = cdf0(ti2,ti1+1) + cdf0(ti2,ti1)
        end do
      end do
      do ti1=1,(Ntrl-1)
        do ti2=1,(Ntrl)
          cdf1(ti2,ti1+1) = cdf1(ti2,ti1+1) + cdf1(ti2,ti1)
        end do
      end do
      ksdiff(4) = maxval(abs(cdf0-cdf1))

      KS(xi) = maxval(ksdiff)
    end do
    !$omp end parallel do

end subroutine
