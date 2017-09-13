subroutine pinput
  use parameters
  use inputFiles
  implicit none
  integer, external :: getpid
    


  write(tmpfile,"(Z5.5)") getpid()
  tmpfile='tmpworkingfile_for_DSM3D'//tmpfile
 
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(5,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  
  open(unit=1,file=tmpfile,status='unknown')
  read(1,110) DSMconfFile
  read(1,110) outputDir
  read(1,110) psvmodel
  read(1,110) modelname
  print *, modelname
  outputDir=trim(outputDir)
  psvmodel=trim(psvmodel)
  modelname=trim(modelname)
  read(1,*) tlen
  read(1,*) rmin_,rmax_,rdelta_
  read(1,*) r0min
  r0max=r0min
  r0delta=20.d0
  read(1,*) thetamin,thetamax,thetadelta
  read(1,*) imin,imax
  read(1,*) rsgtswitch,tsgtswitch,synnswitch
  close(1,status='delete')
 
  ! making directories

  commandline = 'mkdir -p '//trim(outputDir)
  call system(commandline)
  commandline = 'mkdir -p '//trim(outputDir)//'/RSGT'
  call system(commandline)
  commandline = 'mkdir -p '//trim(outputDir)//'/TSGT'
  call system(commandline)
  commandline = 'mkdir -p '//trim(outputDir)//'/log'
  call system(commandline)  


end subroutine pinput
