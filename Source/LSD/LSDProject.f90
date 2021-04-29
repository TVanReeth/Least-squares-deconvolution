PROGRAM LSDProject
  USE FortranGeneral
  USE DataIO
  USE PhysicalConstants, only : c
  USE SpectralAnalysis
  
  IMPLICIT NONE
  
  character :: varfile*(150), resultfile*(150), file100*(150), file101*(150), file102*(150), file103*(150)
  type(dataset) :: spec
  type(dataset), allocatable :: model(:), mod_red(:), ccf(:), lsd(:)
  double precision :: wb, we, regular, binheight, binwidth
  double precision, allocatable :: llim(:), vb(:), ve(:), dv(:), vb2(:), ve2(:), dv2(:), temp(:,:), intg(:)
  integer :: i, j, k, x, t1, t2, clock_rate, clock_max, top, bottom, row, Open_Status, nsec, min, nmin, nh, msec, nlsd, lsdcnt
  
  call System_clock (t1, clock_rate, clock_max)
  
  ! Read in the data from the file UserInput.dat
  call ReadUserInput(spec, model, varfile, resultfile, wb, we, regular, nlsd, llim, vb, ve)
  
  
  ! Read only a part of the spectrum
  call ReadData(trim(spec%name), temp)
  bottom = 1+count(temp(:,1) < wb)
  top = count(temp(:,1) < we)
  allocate(spec%val(top+1- bottom,2))
  spec%val(:,1) = log(temp(bottom:top,1))
  spec%val(:,2) = temp(bottom:top,2)
  deallocate(temp)
  
  spec%n = size(spec%val,1)
  
  ! Read the variance, if provided
  if(trim(varfile)=='./Data/Spectra/') then
      allocate(spec%err(size(spec%val,1)))
      spec%err = 1d0
  else
      call ReadData(trim(varfile), temp)
      allocate(spec%err(size(temp,1)))
      spec%err = sqrt(temp(:,2))
  end if
  
  ! Define the velocity step size dv
  allocate(dv(size(model)))
  dv = c * sum(spec%val(2:,1)-spec%val(:spec%n-1,1)) / (spec%n-1)
  
  ! Read the models
  do i = 1,size(model)
      call ReadData(trim(model(i)%name), model(i)%val)
      model(i)%n = size(model(i)%val,1)
  end do
  
  ! And adapt them to our needs
  call SelectLines(model, mod_red, wb, we, vb, ve, nlsd, llim) 
  
  ! Obtain results
  call LSDCalc(spec, lsd, mod_red, vb, ve, dv, nlsd, regular)
  
  ! Save the obtained results
  file100 = trim(resultfile)//'LSC_noblocks_input.conf'
  file101 = trim(resultfile)//'line_list.asc'
  file102 = trim(resultfile)//'LSD_profile.asc'
  
  
  open(unit=100, file=trim(adjustl(file100)), status="replace", action = "write", IOSTAT = Open_Status)
  
  if (Open_Status /= 0) THEN
      print *,'STATUS: ',Open_Status
      STOP "ERROR: Saving the results: the file `LSC_noblocks_input.conf` was not opened properly."
  end if
  
  write(100,"(A,2x,A)") trim(adjustl(spec%name)),'! Observed spectrum'
  write(100,"(A)") '10 ! Number of iterations'
  write(100,"(A,2x,A)") trim(adjustl(NumStr(size(model)))), '! Number of stellar components'
  write(100,"(A,2x,A)") trim(adjustl(file101)), '! Line lists extension'
  write(100,"(A,2x,A)") trim(adjustl(file102)), '! LSD profiles extension'
  write(100,"(i2,2x,A)") nlsd,'! Number of LSD components'
  do i = 1,nlsd-1
      write(100,"(f5.3,2x)",advance="no") llim(i)
  end do
  write(100,"(A)") '! Theoretical line strength limits (nLSD-1)'
  write(100,"(i5,2x,i5,2x,A)") int(wb),int(we),'! Wavelength range'
  
  close(unit=100)
  
  open(unit=101, file=trim(adjustl(file101)), status="replace", action = "write", IOSTAT = Open_Status)
  if (Open_Status /= 0) then
      print *,'STATUS: ',Open_Status
      STOP "ERROR: Saving the results: the file `line_list.asc` was not opened properly."
  end if
  
  open(unit=102, file=trim(adjustl(file102)), status="replace", action = "write", IOSTAT = Open_Status)
  if (Open_Status /= 0) then
      print *,'STATUS: ',Open_Status
      STOP "ERROR: Saving the results: the file `LSD_profile.asc` was not opened properly."
  end if
  
  do i = 1,size(model)
      file103 = trim(adjustl(resultfile))//'comp'//trim(adjustl(NumStr(i)))//'.lsd'
      open(unit=103, file=trim(adjustl(file103)), status="replace", action = "write", IOSTAT = Open_Status)
      
      if (Open_Status /= 0) then
          print *,'STATUS: ',Open_Status
          STOP "ERROR: Saving the results: the file was not opened properly."
      end if
      
      do row = 1, lsd(nlsd*i)%n-1
          write(103, "(f15.8,2x)",advance="no") lsd(nlsd*i)%val(row,1)
          do j = 1,nlsd
              lsdcnt = nlsd*(i-1)+j
              write(103, "(2x,f12.8)",advance="no") 1.d0-lsd(lsdcnt)%val(row,2) 
          end do
          do j = 1,nlsd-1
              lsdcnt = nlsd*(i-1)+j
              write(103, "(2x,f12.8)",advance="no") lsd(lsdcnt)%err(row)
          end do
          lsdcnt = nlsd*i
          write(103, "(2x,f12.8)") lsd(lsdcnt)%err(row)
      end do
      
      close(unit=103)
      
      write(101,"(A)") trim(adjustl(model(i)%name))
      write(102,"(A)") trim(adjustl(file103))
  end do
  
  close(unit=101)
  close(unit=102)
  
  
  ! Saving the integrals of the different LSD profiles
  open(unit=200, file=trim(resultfile)//'log', status="unknown", action = "write", position='append', IOSTAT = Open_Status)
  write(200,"(A)") '! Integrals of the LSD profiles for the different models and their sum'
  
  allocate(intg(nlsd+1))
  
  do i = 1,size(model)
      intg = 0.d0
      do j = 1,nlsd
          do k = 1,lsd(nlsd*(i-1)+j)%n-1
              binheight = 0.5* ( lsd(nlsd*(i-1)+j)%val(1+k,2) + lsd(nlsd*(i-1)+j)%val(k,2) )
              binwidth = lsd(nlsd*(i-1)+j)%val(1+k,1) - lsd(nlsd*(i-1)+j)%val(k,1)
              intg(j) = intg(j) + binheight * binwidth
          end do
      end do
      intg(nlsd+1) = sum(intg(1:nlsd))
      write(200,"(10(f15.9,2x))") (intg(j), j=1,nlsd+1)
  end do
  
  close(unit=200)
  
  deallocate(model,mod_red,lsd, intg)
  
  call System_clock(t2, clock_rate, clock_max)
  
  nsec = int(real( t2 - t1 ) / real(clock_rate))
  nmin = nsec/60; nh = nmin/60
  min  = nmin - nh*60
  msec = nsec - min*60 - nh*3600
  write(*,*)
  write(*,'(a,i8,a,i8,a,i8,a)') ' Calculation time: ', nh, ' h', min, ' min', msec, ' sec'
  
END PROGRAM LSDProject
