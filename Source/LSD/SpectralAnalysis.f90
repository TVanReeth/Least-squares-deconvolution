! -*- f90 -*-
MODULE SpectralAnalysis
!==============================================================================

  IMPLICIT NONE

  
  interface LSDCalc
      subroutine LSDCalc(spec, lsd, model, vb, ve, dv, nlsd, regular)
          use FortranGeneral
          implicit none
          type(dataset), intent(in) :: spec, model(:)
          type(dataset), allocatable, intent(inout) :: lsd(:)
          double precision, intent(in) :: vb(:), ve(:), dv(:), regular
          integer, intent(in) :: nlsd
      end subroutine LSDCalc
  end interface
  
  
  interface SelectLines
      subroutine SelectLines(model, mod_red, wb, we, vb, ve, nlsd, llim)
          USE FortranGeneral
          implicit none
          type(dataset), intent(in) :: model(:)
          type(dataset), allocatable, intent(inout) :: mod_red(:)
          double precision, intent(in) :: wb, we, vb(:), ve(:), llim(:)
          integer, intent(in) :: nlsd
      end subroutine SelectLines
  end interface
  
  
  interface SynthSpec
      subroutine SynthSpec(modlines, lsdprof, synth, n)
          use FortranGeneral
          implicit none
          type(dataset), intent(in) :: modlines(:), lsdprof(:)
          integer, intent(in) :: n
          double precision, intent(inout) :: synth(n,2)
      end subroutine SynthSpec
  end interface

END MODULE SpectralAnalysis
 
 
 
SUBROUTINE LSDCalc(spec, lsd, model, vb, ve, dv, nlsd, regular)
     !Calculate the LSD profile
     USE FortranGeneral
     USE MathRoutines
     USE PhysicalConstants, only : c
     USE DataIO
     IMPLICIT NONE
     
     type(dataset), intent(in) :: spec, model(:)
     double precision, intent(in) :: vb(:), ve(:), dv(:), regular
     integer, intent(in) :: nlsd
     type(dataset), allocatable, intent(inout) :: lsd(:)
     integer :: i,j,k,l, tot, wb, we, pt, tot2
     double precision :: f1, meanmod, left(2), right(2), center(2)
     double precision, allocatable :: M(:,:), S2(:), spec_red(:,:), tmp(:), tmp_err(:), velocity(:)
     
     
     allocate(lsd(size(model)))
     
     tot = 0
     
     do l=1,size(model)
         lsd(l)%n = int(1 + (ve(ceiling(dble(l)/dble(nlsd))) - vb(ceiling(dble(l)/dble(nlsd)))) / dv(ceiling(dble(l)/dble(nlsd))))
         allocate(lsd(l)%val(lsd(l)%n,2))
         allocate(lsd(l)%err(lsd(l)%n))
         lsd(l)%val(:,1) = (/(vb(ceiling(dble(l)/dble(nlsd))) + (j-1)*dv(ceiling(dble(l)/dble(nlsd))), j=1,lsd(l)%n)/)
         tot = tot + lsd(l)%n + 2 ! Include a ghost pixel on either side of each LSD
         
     end do
     
     allocate(spec_red(spec%n,2))
     spec_red = spec%val
     
     allocate(S2(spec%n))
     !S2 = 1d0/(spec%n**2)
     S2 = spec%n / sum((spec%val(:,2)-sum(spec%val(:,2))/spec%n)**2.d0)
     
     do l = 1,spec%n
         if(spec%err(l) /= 0) then
             S2(l) = 1d0/(spec%err(l)**2)
         else
             spec_red(l,2) = 0d0
         end if
     end do
     
     S2 = S2/sum(S2)
     
     spec_red(:,2) = spec_red(:,2) - sum(S2 * spec_red(:,2))
     
     do l = 1,spec%n
         if(spec%val(l,2) == 0d0) then
             spec_red(l,2) = 0d0
         end if
     end do
     
     allocate(M(spec%n,tot))
     M = 0
     
     pt = 0
     
     meanmod = 0
     
     do l=1,size(model)
         
         if(nlsd > 1) then
             if(mod(l,nlsd) == 1) then
                 meanmod = 0.d0
                 tot2    = 0.d0
                 do k = 1,nlsd
                     meanmod = meanmod + sum(model(l+k-1)%val(:,2))
                     tot2    = tot2 + model(l+k-1)%n
                 end do
                 meanmod = meanmod / tot2
             end if
         else
             meanmod = sum(model(l)%val(:,2)) / model(l)%n
         end if
         
         do k=1,model(l)%n
             allocate(velocity(lsd(l)%n+2))
             velocity(1) = 2*lsd(l)%val(1,1) - lsd(l)%val(2,1)
             velocity(2:lsd(l)%n+1) = lsd(l)%val(:,1)
             velocity(lsd(l)%n+2) = 2*lsd(l)%val(lsd(l)%n,1) - lsd(l)%val(lsd(l)%n-1,1)
             wb = 1+count(spec_red(:,1) < model(l)%val(k,1) + velocity(1)/c)
             we = count(spec_red(:,1) < model(l)%val(k,1) + velocity(lsd(l)%n+2)/c)
             do i=wb,we
                 j = count(velocity .le. c*(spec_red(i,1)-model(l)%val(k,1)))
                 f1 = (c*(spec_red(i,1)-model(l)%val(k,1)) - velocity(j))/(velocity(j+1)-velocity(j))
                 M(i,pt+j) = M(i,pt+j) + (model(l)%val(k,2)-meanmod)*(1-f1)
                 M(i,pt+j+1) = M(i,pt+j+1) + (model(l)%val(k,2)-meanmod)*f1
             end do
             deallocate(velocity)
         end do
         
         pt = pt + lsd(l)%n + 2
         
     end do
     
     allocate(tmp(tot), tmp_err(tot))
     
     call RegLS(M, S2, spec_red(:,2), tmp, tmp_err, regular) ! Solve regularized Y0 = M.x, with weights S2, using LSD
     
     deallocate(M, S2)
     
     pt = 0
     
     do l=1,size(lsd)
         lsd(l)%val(:,2) = tmp(pt+2:pt+lsd(l)%n)
         lsd(l)%err = tmp_err(pt+2:pt+lsd(l)%n)
         pt = pt + lsd(l)%n + 2 ! Take the ghost pixels on either side of each LSD profile into account !
     end do
     
     deallocate(tmp,spec_red)
     
END SUBROUTINE LSDCalc



SUBROUTINE SelectLines(model, mod_red, wb, we, vb, ve, nlsd, llim)
    USE MathRoutines
    USE FortranGeneral
    USE DataIO
    USE PhysicalConstants
    implicit none
    type(dataset), intent(in) :: model(:)
    type(dataset), allocatable, intent(inout) :: mod_red(:)
    double precision, intent(in) :: wb, we, vb(:), ve(:), llim(:)
    integer, intent(in) :: nlsd
    integer :: nr, l, i, j, k
    double precision :: blim, tlim
    type(dataset), allocatable :: tmp_red(:)
    
    
    nr = size(model)
    
    allocate(tmp_red(nr))
    allocate(mod_red(nlsd*nr))
    
    do l=1,nr
            
        ! Make sure the model only includes lines which fall within the borders of 
        ! the given spectrum
        tmp_red(l)%n = count( (log(wb)-vb(l)/c .le. log(model(l)%val(:,1))) .and. (log(model(l)%val(:,1)) .le. log(we)-ve(l)/c) )
        
        allocate(tmp_red(l)%val(tmp_red(l)%n,2))
        
        k = 1
        
        do i = 1,model(l)%n
            if((log(wb)-vb(l)/c .le. log(model(l)%val(i,1))) .and. (log(model(l)%val(i,1)) .le. log(we)-ve(l)/c)) then
                tmp_red(l)%val(k,1) = log(model(l)%val(i,1))
                tmp_red(l)%val(k,2) = model(l)%val(i,2)
                k = k+1
            end if
        end do
        
        ! Split the theoretical lines in the models in different groups according to the line strengths
        
        if(nlsd > 1) then
            do i = 1, nlsd
                
                mod_red(nlsd*(l-1)+i)%name = trim(model(l)%name)//'_'//trim(adjustl(NumStr(i)))
                
                if(i == 1) then
                    tlim = 1.d0
                else
                    tlim = llim(nlsd+1-i)
                end if
                
                if(i == nlsd) then
                    blim = 0.d0
                else
                    blim = llim(nlsd-i)
                end if
                
                mod_red(nlsd*(l-1)+i)%n = count((tmp_red(l)%val(:,2) .ge. blim) .and. (tmp_red(l)%val(:,2) < tlim))
                
                allocate(mod_red(nlsd*(l-1)+i)%val(mod_red(nlsd*(l-1)+i)%n,2))
                k = 1
                
                do j = 1,tmp_red(l)%n
                    if((tmp_red(l)%val(j,2) .ge. blim) .and. (tmp_red(l)%val(j,2) < tlim)) then
                        mod_red(nlsd*(l-1)+i)%val(k,:) = tmp_red(l)%val(j,:)
                        k = k+1
                    end if
                end do
                
            end do
        
        else
            
            mod_red(l)%name = model(l)%name
            mod_red(l)%n = size(tmp_red(l)%val,1)
            allocate(mod_red(l)%val(mod_red(l)%n,2))
            mod_red(l)%val = tmp_red(l)%val
            
        end if
        
    end do
    
    deallocate(tmp_red)

END SUBROUTINE SelectLines



SUBROUTINE SynthSpec(modlines, lsdprof, synth, n)
    USE FortranGeneral
    USE MathRoutines
    USE PhysicalConstants
    IMPLICIT NONE
    type(dataset), intent(in) :: modlines(:), lsdprof(:)
    integer, intent(in) :: n
    double precision, intent(inout) :: synth(n,2)
    integer :: i, j
    double precision :: line, depth!, Interp(n)
    
    synth(:,2) = 0.
    
    do i = 1,size(lsdprof)
        do j = 1,modlines(i)%n
            line = modlines(i)%val(j,1)
            depth = 1 - modlines(i)%val(j,2)
            synth(:,2) = synth(:,2) + depth * Interp(log(synth(:,1)), line + lsdprof(i)%val(:,1)/c, lsdprof(i)%val(:,2))
        end do
    end do
    
    synth(:,2) = 1d0 - synth(:,2)
    
END SUBROUTINE SynthSpec
