MODULE FortranGeneral
  IMPLICIT NONE
  PRIVATE
  
  public NumStr, dataset, indexlist
  !! Various handy python-like functions, which are to advanced to be included in standard Fortran 90/95.
  
  ! Variables dataset and indexlist, created to facilitate listing the different models. (I.e. equivalent of a special case of python's "list")
  type dataset
      character(len=150) :: name
      integer :: n = 0
      double precision, allocatable :: val(:,:), err(:)
  end type
  
  type indexlist
      integer, allocatable :: val(:)
  end type
  
  ! Converting numbers of different types to a string
  interface NumStr
      function IntStr(number)
          implicit none
          integer :: number
          character :: IntStr*(24)
      end function
      
      function RealStr(number)
          implicit none
          real :: number
          character :: RealStr*(24)
      end function
      
      function DbleStr(number)
          implicit none
          double precision :: number
          character :: DbleStr*(24)
      end function
  end interface

END MODULE


function IntStr(number)
    implicit none
    integer :: number
    character :: IntStr*(24) 
    
    write(IntStr,*) number
    
end function


function RealStr(number)
    implicit none
    real :: number 
    character :: RealStr*(24)
    
    write(RealStr,*) number
    
end function


function DbleStr(number)
    implicit none
    double precision :: number
    character :: DbleStr*(24)
    
    write(DbleStr,*) number
    
end function