! -*- f90 -*-
MODULE MathRoutines
!==============================================================================
  IMPLICIT NONE
  !! Various handy mathematical routines
  
  ! Derive the (function) values in an array, which are mapped onto a certain domain.
  interface Derive
      subroutine Derive(domain, array, deriv)
          implicit none
          double precision, intent(in) :: domain(:), array(:)
          double precision, allocatable, intent(inout) :: deriv(:)
      end subroutine Derive
  end interface
  
  ! Compute the inverse of a matrix, using LU-decomposition
  interface Inverse
      subroutine Inverse(matrix, inv_mat)
          implicit none
          double precision, dimension(:,:), intent(in) :: matrix
          double precision, allocatable, intent(inout) :: inv_mat(:,:)
      end subroutine Inverse
  end interface
  
  ! Compute the regularised solution of the matrix equation Y0 = M*x, with S2 the variance of the observed solution Y0, i.e.
  ! x = (M*S2*M)^(-1) * (M*S2*Y0). Here regular is the regularization parameter, 
  !and R is the regularization matrix. (See also Kochukhov et al. (2010))
  interface RegLS
      subroutine RegLS(M, S2, Y0, x, x_err, regular)
          implicit none
          double precision, intent(in) :: M(:,:), S2(:), Y0(:), regular
          double precision, intent(inout) :: x(:), x_err(:)
      end subroutine RegLS
  end interface
  
  ! Give the index of the first element which is true in a boolean array
  interface Firstloc
      function FirstLoc(array)
          implicit none
          logical, dimension(:) :: array
          integer :: FirstLoc
      end function FirstLoc
  end interface
  
  ! Spline & splint: adapted from numerical recipies
  interface Spline
      subroutine Spline(xn,yn, y2n)
          implicit none
          double precision, dimension(:), intent(in) :: xn, yn
          double precision, dimension(:),allocatable, intent(inout) :: y2n
      end subroutine Spline
  end interface
  
  interface Splint
      subroutine Splint(x,y,xn,yn,y2n)
          implicit none
          double precision, dimension(:), intent(in) :: x,xn,yn,y2n
          double precision, dimension(:), intent(inout) :: y
      end subroutine Splint
  end interface
  
  ! Sorting algorithm
  interface SSort
      subroutine SSort(x, y)
          implicit none
          double precision, intent(inout) :: x(:), y(:)
      end subroutine SSort
  end interface
  
  
  contains
  
  ! Linear interpolation, similar to numpy.interp() in python
  function Interp(x, xp, fp, left, right)
        IMPLICIT NONE
        double precision, intent(in) :: x(:), xp(:), fp(:)
        double precision, optional, intent(in) :: left, right
        double precision :: xb, xe, f1
        double precision,allocatable :: Interp(:)
        integer :: i, j, k, l
        
        if(present(left)) then
            xb = left
        else
            xb = 0
        end if
        
        if(present(right)) then
            xe = right
        else
            xe = 0
        end if
        
        allocate(Interp(size(x)))
        
        Interp = 0
        
        i = count(x < xp(1))
        j = count(x .le. xp(size(xp)))
        
        if(i > 0) then
            Interp(:i) = xb
        end if
        
        if(j < size(x)) then
            Interp(j+1:) = xe
        end if
        
        do k = i+1, j
            l = count(xp < x(k))
            f1 = (x(k) - xp(l))/(xp(l+1) - xp(l))
            Interp(k) = fp(l)*(1-f1) + fp(l+1)*f1
        end do
        
        return
        
        deallocate(Interp)
        
  end function Interp
  
  
  ! Quadratic interpolation function, provided to me by Andrew Tkachenko
  function map1 (xold, fold, nold, xnew, fnew, nnew)
      implicit none
      double precision :: xold(*), fold(*), xnew(*), fnew(*)
      integer :: l, l1, l2, ll, k, nnew, nold, map1
      double precision :: a, b, c, d, abac, afor, bbac, bfor, cbac, cfor, wt
      
      l  = 2
      ll = 0
      do k = 1, nnew
          do while ( xnew(k) >= xold(l) .and. l <= nold )
              l = l + 1
          enddo
          if ( l > nold ) then
              if ( l /= ll ) then
                  
                  l  = amin0(nold,l)
                  c  = 0.d0
                  b  = (fold(l) - fold(l-1)) / (xold(l) - xold(l-1))
                  a  = fold(l) - xold(l) * b
                  ll = l
                  
              end if
              
              fnew(k) = a + (b + c * xnew(k)) * xnew(k)
              cycle
          end if
          if ( l /= ll ) then
              if ( l /= 2 .and. l /= 3 ) then
                  l1 = l - 1
                  if ( .not. ( l > ll+1 .or. ( l == 3 .or. l == 4 ) ) ) then
                      
                      cbac = cfor
                      bbac = bfor
                      abac = afor
                      if ( l == nold ) then
                          
                          c  = cbac
                          b  = bbac
                          a  = abac
                          ll = l
                          fnew(k) = a + (b + c * xnew(k)) * xnew(k)
                          cycle
                          
                      end if
                  else
                      
                      l2 = l - 2
                      d  = (fold(l1) - fold(l2)) / (xold(l1) - xold(l2))
                      cbac = fold(l) / ((xold(l) - xold(l1)) * (xold(l) - xold(l2))) + &
                             (fold(l2) / (xold(l) - xold(l2)) - fold(l1) /             &
                             (xold(l) - xold(l1))) / (xold(l1) - xold(l2))
                      bbac = d - (xold(l1) + xold(l2)) * cbac
                      abac = fold(l2) - xold(l2) * d + xold(l1) * xold(l2) * cbac
                      if ( l >= nold ) then
                          
                          c  = cbac
                          b  = bbac
                          a  = abac
                          ll = l
                          fnew(k) = a + (b + c * xnew(k)) * xnew(k)
                          cycle
                          
                      endif
                  endif
                  d = (fold(l) - fold(l1)) / (xold(l) - xold(l1))
                  cfor = fold(l+1) / ((xold(l+1) - xold(l)) * (xold(l+1) - xold(l1))) + &
                         (fold(l1) / (xold(l+1) - xold(l1)) - fold(l) /                 &
                         (xold(l+1) - xold(l))) / (xold(l) - xold(l1))
                  bfor = d - (xold(l) + xold(l1)) * cfor
                  afor = fold(l1) - xold(l1) * d + xold(l) * xold(l1) * cfor
                  wt = 0.d0
                  if ( abs(cfor) /= 0.d0 ) wt = abs(cfor) / ( abs(cfor) + abs(cbac) )
                  
                  a  = afor + wt * (abac - afor)
                  b  = bfor + wt * (bbac - bfor)
                  c  = cfor + wt * (cbac - cfor)
                  ll = l
                  fnew(k) = a + (b + c * xnew(k)) * xnew(k)
                  cycle
              endif
              if ( l /= ll ) then
                  
                  l  = amin0(nold,l)
                  c  = 0.d0
                  b  = (fold(l) - fold(l-1)) / (xold(l) - xold(l-1))
                  a  = fold(l) - xold(l) * b
                  ll = l
                  
              endif
          endif
          fnew(k) = a + (b + c * xnew(k)) * xnew(k)
      enddo
      map1    = ll - 1
  end function map1
  
END MODULE MathRoutines


subroutine Derive(domain, array, deriv)
    implicit none
    double precision, intent(in) :: domain(:), array(:)
    double precision, allocatable, intent(inout) :: deriv(:)
    
    allocate(deriv(size(array)))
    deriv(2:size(array)-1) = (array(3:)-array(:size(array)-2))/(domain(3:) - domain(:size(domain)-2))
    deriv(1) = (array(2) - array(1))/(domain(2) - domain(1))
    deriv(size(array)) = (array(size(array)) - array(size(array)-1))/(domain(size(domain)) - domain(size(domain)-1))
end subroutine Derive


subroutine Inverse(matrix, inv_mat)
    implicit none
    double precision, dimension(:,:), intent(in) :: matrix
    double precision, allocatable, intent(inout) :: inv_mat(:,:)
    integer :: i,j,n, nb, lwork, info, ilaenv
    integer, allocatable :: ipiv(:)
    double precision :: cond
    double precision, allocatable :: work(:)
    external :: dgetrf, ilaenv, dgetri
    
    if(size(matrix,1) /= size(matrix,2)) then
        STOP "ERROR: MathRoutines: Inverse: The given matrix is not squared."
    end if
    
    n = size(matrix,1)
    allocate(inv_mat(n,n))
    allocate(ipiv(n))
    inv_mat = matrix
    call dgetrf(n, n, inv_mat, n, ipiv, info)
    
    if(info /= 0) then
        print *, "info = ",info
        STOP "ERROR: MathRoutines: Inverse: occurred during the LU-decomposition."
    end if
    
    nb = ilaenv(1, 'MathRoutines', ' ', n, -1, -1, -1)
    
    if(nb < 0) then
        print *,"NB = ",nb
        STOP "ERROR: MathRoutines: Inverse: Wrong argument for the function ILAENV."
    end if
    
    lwork = n*nb
    allocate(work(lwork))
    
    call dgetri(n, inv_mat, n, ipiv, work, lwork, info)
    
    if(info /= 0)  STOP "ERROR: MathRoutines: Inverse: occurred during the matrix inversion."
    
    deallocate(work)
    deallocate(ipiv)
    
    cond = maxval(sum(matrix,1))*maxval(sum(inv_mat,1))
    print *, "MathRoutines: Inverse: condition number = ",cond

end subroutine Inverse



subroutine RegLS(M, S2, Y0, x, x_err, regular)
    USE DataIO
    USE MathRoutines, only : Inverse
    implicit none
    double precision, intent(in) :: M(:,:), S2(:), Y0(:), regular
    double precision, intent(inout) :: x(:), x_err(:)
    double precision, allocatable :: MATR(:,:), cov(:,:), tmp(:), MS2T(:,:), A(:,:), b(:)
    integer :: i,j,l,n 
    double precision :: lb, le, step, resc
    
    
    n = size(x)
    
    ! Construct the matrix R
    allocate(MATR(n,n))
    
    MATR = 0
    
    do i=1,n
        if(i < n) then
            MATR(i,i+1) = -1
            MATR(i+1,i) = -1
        end if
        if (i < n-1) MATR(i+1,i+1) = 2
    end do
    
    MATR(1,1) = 1
    MATR(n,n) = 1
    
    allocate(MS2T(n,size(M,1)))
    
    do i=1,size(M,1)
        MS2T(:,i) = M(i,:)*S2(i)
    end do
    
    allocate(A(n,n))
    allocate(b(n))
    
    A = matmul(MS2T,M)
    b = matmul(MS2T,Y0)
    
    deallocate(MS2T)
    
    allocate(tmp(n))
    
    call Inverse(A,cov)
    tmp = b*maxval(matmul(cov,b))/maxval(b)
    deallocate(cov)
    
    call Inverse(A + regular*MATR, cov)
    x = matmul(cov, b + matmul(regular*MATR,tmp))
!     x = matmul(cov, b)
    
    deallocate(MATR, A, b, tmp)
    
    !TODO: Check the error calculations !!!
    do i=1,n
        x_err(i) = sqrt(cov(i,i) * sum(S2*(Y0 - matmul(M,x))**2)/(size(Y0)-n))
    end do
    
    deallocate(cov)
    
end subroutine RegLS


function FirstLoc(array)
    implicit none
    logical, dimension(:) :: array
    integer :: FirstLoc, i
    
    FirstLoc = 0
    
    do i=1,size(array)
        if(array(i)) then
            FirstLoc = i
            exit
        end if
    end do
    
    return
    
end function Firstloc



subroutine Spline(xn,yn,y2n)
    ! adapted from "Numerical Recipes in Fortran 77, 2nd ed." 
    ! by W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery
    implicit none
    double precision, dimension(:), intent(in) :: xn, yn
    double precision, dimension(:), allocatable, intent(inout) :: y2n
    double precision :: yp1, ypn, p, qn, sig, un
    double precision, allocatable :: u(:)
    integer :: i,n
    
    n = size(xn)
    allocate(u(n))
    allocate(y2n(n))
    
    yp1 = (yn(2)-yn(1)) / (xn(2) - xn(1))
    ypn = (yn(n)-yn(n-1)) / (xn(n) - xn(n-1))
    
    if(yp1 .gt. .99e30) then
        y2n(1) = 0.
        u(1) = 0.
    else
        y2n(1) = -0.5
        u(1) = (3./(xn(2)-xn(1)))*((yn(2)-yn(1))/(xn(2)-xn(1))-yp1)
    end if
    
    do i=2,n-1
        sig = (xn(i)-xn(i-1)) / (xn(i+1)-xn(i-1))
        p = sig*y2n(i-1) + 2.
        y2n(i) = (sig-1.)/p
        u(i) = ( 6. * ((yn(i+1)-yn(i))/(xn(i+1)-xn(i)) - (yn(i)-yn(i-1))/(xn(i)-xn(i-1))) / (xn(i+1)-xn(i-1)) - sig*u(i-1))/p
    end do
    
    if(ypn .gt. .99e30) then
        qn = 0.
        un = 0.
    else
        qn = 0.5
        un = (3./(xn(n)-xn(n-1))) * (ypn-(yn(n)-yn(n-1))/(xn(n)-xn(n-1)))
    end if
    
    y2n(n) = (un - qn*u(n-1)) / (qn*y2n(n-1)+1.)
    
    do i = n-1, 1, -1
        y2n(i) = y2n(i) * y2n(i+1) + u(i)
    end do
    
    deallocate(u)
    
    return
    
end subroutine Spline



subroutine Splint(x, y, xn, yn, y2n)
    ! taken from "Numerical Recipes in Fortran 77, 2nd ed." 
    ! by W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery
    implicit none
    double precision, dimension(:), intent(in) :: x, xn, yn, y2n
    double precision, dimension(:), intent(inout) :: y
    integer :: i, k, klo, khi
    double precision :: h, a, b
    
    
    do i = 1, size(x)
        
        klo = 1
        khi = size(xn)
        
    1   if(khi-klo .gt. 1) then
            k = (khi+klo)/2.
            if(xn(k) .gt. x(i)) then
                khi = k
            else
                klo = k
            end if
        goto 1
        end if
        
        h = xn(khi)-xn(klo)
        if(h .eq. 0.) pause 'MathRoutines: Splint: bad xn input in splint'
        
        a = (xn(khi)-x(i))/h
        b = (x(i)-xn(klo))/h
        
        y(i) = a*yn(klo) + b*yn(khi) + ((a**3-a)*y2n(klo)+(b**3-b)*y2n(khi))*(h**2)/6.
        
    end do
    
    return
    
end subroutine Splint


SUBROUTINE SSort(x, y)
    IMPLICIT NONE
    double precision, intent(inout) :: x(:), y(:)
    double precision :: temp
    integer :: n, i, j, k, itemp
    
    
    n = size(x)
    
    do i=2,n
        
        if (x(i).lt.x(i-1)) then
            do j=i-2,1,-1
                    if(x(i).gt.x(j)) go to 70
            end do
            j=0
            
  70        x(j+1:i) = cshift(x(j+1:i),-1)
            y(j+1:i) = cshift(y(j+1:i),-1)
            
        end if
    end do
    return
END SUBROUTINE SSort