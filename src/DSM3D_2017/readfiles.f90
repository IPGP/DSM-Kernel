subroutine parameterInput
  use DSMparamters
  use parameters
  use inputFiles
  implicit none
  integer, external :: getpid
  character(200) :: firstinf
  character(120) :: dummy
  character(120) :: commandline
  character(200) :: model1d


  call getarg(1,firstinf)
  firstinf=trim(firstinf)

  write(jobid,"(Z5.5)") getpid()
  tmpfile='tmpworkingfile_for_DSM3D'//jobid
 
  open(unit=5, file=firstinf,status='unknown')
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
  close(5)
  
  open(unit=1,file=tmpfile,status='unknown')
  read(1,110) DSMconfFile
  read(1,110) outputDir
  read(1,110) psvmodel
  read(1,110) modelname
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


  call readDSMconf
  if(index(psvmodel,".card")) then ! "psvmodel" should be interpreted as a normal mode cardif ".card" is in its string
     model1d=trim(psvmodel)


  endif



end subroutine parameterInput




subroutine readDSMconf
  use inputFiles
  use DSMparameters
  implicit none
  character(120) :: dummy
  character(120) :: tmpfile

  tmpfile='tmpworkingfile_for_DSMconf'//jobid

  open(unit=2, file=DSMconfFile, status='old',action='read',position='rewind')
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(2,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(2)
  
 
  open(unit=1,file=tmpfile,status='unknown')
  read(1,*) re
  read(1,*) ratc
  read(1,*) ratl
  read(1,*) omegai
  read(1,*) maxlmax
  close(1,status='delete')


end subroutine readDSMconf
  

subroutine readpsvmodel(psvmodel,tmpfile)
  implicit none
  character(120) :: psvmodel, tmpfile, dummy
  open(unit=2, file=psvmodel, status='old',action='read',position='rewind')
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(2,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(2)
end subroutine readpsvmodel
