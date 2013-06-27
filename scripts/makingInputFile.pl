#!/usr/bin/perl


$DSMHOME = "/home/dtp02/fuji/DSM-KS2011/";
$modelname = "ak135";
$outputDirDir = "/export/bacchus1/fuji/".$modelname;
$npmpi = 8;
unless(-d $outputDirDir){
    system "mkdir $outputDirDir";
}


#--------------------------------------------------------------

$psvmodel= $DSMHOME."/etc/models/".$modelname.".psv.model";
$mpiSGTpsv = $DSMHOME."/bin/mpiSGTpsv";
$mpiSGTsh  = $DSMHOME."/bin/mpiSGTsh";

$switching = "1 1 1"; # for input.inf #9
$inf_calculDir = $DSMHOME."/inf_calcul/";
$DSMconfig = $DSMHOME."/etc/DSMconfiguration/DSM.config";

for(my $i=1;$i<=1;$i++){
    
    $depthkm=10.e0*$i;
    
    

#--------------------------------------------------------------
    
    $depth=sprintf '%d',$depthkm;
    
    $outputDir =$outputDirDir."/depth".$depth;
    $logDir = $outputDir."/log";
    $RSGTDir = $outputDir."/RSGT";
    $TSGTDir = $outputDir."/TSGT";
    $inputFile = $inf_calculDir."/input".$depth.".inf";
    $qsubFile =  $inf_calculDir."/qsub.".$depth;
    $execFile = $inf_calculDir."/exec".$depth.".sh";
    $logFile = $DSMHOME."/log/".$depth.".".$modelname.".psv.log";
    $logFileSH = $DSMHOME."/log/".$depth.".".$modelname.".sh.log";
    unless(-d $outputDir){
	system "mkdir $outputDir";
    }
    unless(-d $logDir){
	system "mkdir $logDir";
    }
    unless(-d $RSGTDir){
	system "mkdir $RSGTDir";
    }
    unless(-d $TSGTDir){
	system "mkdir $TSGTDir";
    }

    
    open FILE, ">$inputFile";
    print FILE "\# This is the input file for SGTforPinv generated by makingInputFile.pl\n";
    print File "\#\n";
    print FILE "\# 0. DSM configuration file\n";
    print FILE "$DSMconfig \n";
    print FILE "\#\n";
    print FILE "\# 1. outputDir\n";
    print FILE "$outputDir\n";
    print FILE "\# 2. psvmodel\n";
    print FILE "$psvmodel\n";
    print FILE "\# 3. modelname\n";
    print FILE "$modelname\n";
    print FILE "\# 4. tlen (preferred to be set as 2\^n\*0.1 seconds\n";
    print FILE "1.6384d3\n";
    print FILE "\# 5. minQ radius (km), maxQ radius, delta radius\n";
    print FILE "3651.d0 6371.d0 20.d0\n";
    print FILE "\# 6. source radius (km)\n";
    printf FILE '%e', 6371.e0-$depthkm;
    print FILE "\n";
    print FILE "\# 7. minQ theta (degrees), maxQ theta, delta theta\n";
    print FILE "0.d0 180.d0 2.d-1\n";
    print FILE "\# 8. imin imax (omega index: tlen/imax is the shortest period)\n";
    print FILE "0 1024\n";
    print FILE "\# 9. switching calculation of RSGT, TSGT, synthetics\n";
    print FILE "\# if 1 we calculate\n";
    print FILE "$switching\n";
    print FILE "\# don't forget to write \'end\'\n";
    print FILE "end\n";
    close FILE;
    
#exit;
    open FILE, ">$qsubFile";
    print FILE "\#PBS -S /bin/bash\n";
    print FILE "\#PBS -N d$depth\n";
    print FILE "\#PBS -l select=16:ncpus=8:mpiprocs=8\n";
    print FILE "\#PBS -l walltime=100:00:00\n";
    print FILE "\#PBS -j oe\n";
    print FILE "\n";
    print FILE "cd \$PBS_O_WORKDIR \# on se place dans le repertoire courant\n";
    print FILE ". \$MODULESHOME/init/sh\n";
    print FILE "module load default-intel-mpi\n";
    print FILE "cat \$PBS_NODEFILE | uniq \> mpd.hosts$depth\n";
    print FILE "cat \$PBS_NODEFILE \> nodes$depth\n";
    print FILE "mpdboot --rsh=ssh -v -n \`cat mpd.hosts$depth | wc -l\` -f mpd.hosts$depth\n";
    print FILE "mpiexec -machinefile nodes$depth -np $npmpi  $mpiSGTpsv < $inputFile > $logFile\n";
    
    close FILE;
    
    open FILE, ">$execFile";
    print FILE "\#\!/bin/bash\n";
    print FILE "mpirun -np $npmpi $mpiSGTpsv < $inputFile > $logFile\n";
    print FILE "mpirun -np $npmpi $mpiSGTsh  < $inputFile > $logFileSH\n";
    close FILE;
    
    system "chmod +x $execFile";
    
    
    
    #system "qsub $qsubFile";
    

}
