#!/usr/bin/perl

# change models

# Nobuaki Fuji, IPGP, 12/2016
# Nobuaki Fuji, IPGP, 08/2017

# setting local parameters for deconvolution

$param{CUTFREQ1} = 1/360;#STS1
$param{CUTFREQ2} = 1/100;#STS2
$param{CUTFREQDAT} = 1/200;
$param{NYQFREQ} = 10.0;
$param{PI} =  3.1415926535897932;
$param{tapercode} = 1; # 0: no taper, 1: sine, 2: cosine
$param{tw} = 10; #reciprocal of a taper window
$param{DELTA} = 0.05;
$param{deconvolveornot}=1;#if this value equals 0 we don't deconvolve any seismograms

# preventing errors

$param{logKHOLE} = "khole_log.err"; #if this = "000", we do not prevent
$param{logCMPINC} = "cmpinc_log.err"; # if this = "000", we do not prevent
$param{logmerge} = "merge_log.err";
$param{logrotation} = "rotation_log.err";
$param{NEor123}= "NE";
#$param{'logCMPINC'} = "000";

# for eventdata use

$param{zerointerpolating} = 1; # 0: we don't use eventdata, 1: we use them
$param{eventdatatapertime} = 60; # seconds
 
# setting local parameters for DSM calculation

#$param{tlen} = 6553.6;
$param{tlen} = 1638.4;
$param{npts_required} = 0; #if this equals 0, npts will be determined automatically
$param{up2period_syn} = 3.2; #you have to decide, the shortest period you wanna calculate
$param{up2period_par} = 3.2; #for partial derivatives
$param{re} = "1.d-2";
$param{ratc} = "1.d-10";
$param{ratl} = "1.d-5";
$param{wrap} = "1.d-3"; #for the moment we don't calculate separately

# DSM switching

$param{dsmtish} = 0; #don't calculate = 0; calculate = 1
$param{dsmtipsv} = 0;
$param{dsmsshshi} = 0;

# shouldbediscarded

$param{shouldbesiscarded} = 0; # don't delete = 0; delete = 1;

# models

local @model = ("flatshcorrect");
local $modeldir = "../model/";
#ak135 (isotropic), isoPREM (isotropic) and PREM (anisotropic) structure files are not required but otherwise you have to set xxx.sh.model and/or xxx.psv.model file
#but if you put xxx.sh.model and/or xxx.psv.model in the directory $dir, you need not input names of models here! -- searchingModel

local $dir = "./";
local $parentdir = "./";    #if you wanna specify
local $homedir = "./";      #if you wanna specify
local $ijob=0; #never write here



#---------------------------------------------------------------------------


opendir DIR, "$dir" or die "could not open $dir: $!";

local @psvinffiles = grep {/^PSV_/} readdir DIR;

closedir DIR;

opendir DIR, "$dir" or die "could not open $dir: $!";

local @shinffiles = grep {/^SH_/} readdir DIR;

closedir DIR;


mkdir ("../inffiles4malbec", 0755) unless (-e "../inffiles4malbec");
mkdir ("../spcfiles", 0755) unless (-e "../spcfiles");
mkdir ("../trash",0755) unless (-e "../trash");
mkdir ("../submissionfiles",0755) unless(-e "../submissionfiles");



foreach (@psvinffiles){
    $modelname=$_; 
    open FILE, $_;
    local @lines =<FILE>;
    close FILE;
    
    $modelname=~s/\s*$//;

    $modelname=~s/^PSV_//;

    $modelname=~s/\.inf$//;

    #mkdir ("../spcfiles/$_", 0755) unless (-e "../spcfiles/$_");
    #$tishfile="../inffiles/".$_.".sh.inf";
    #$parfile="../inffiles/".$_.".par.inf";

    $tipsvfile="../inffiles4malbec/$_";

    $_=$tipsvfile;

    open FILE1, ">$tipsvfile";

       
    my $lineStart=0;
    my $lineEnd=0;
    my $lineIndex=0;
    my $lineMoments=0;
    foreach(@lines){
	if($lineIndex==17){
	    $_="c";
	}
	$_ =~s/\s*$//;
	#$_ =~s/^..\/spcresult\/model\/1/..\/spcfiles/;
	#$_ =~s/SH.spc$/sh.spc/;
	#print "$_\n";
	if($_=~/.spc$/){
	    my @tmplines=split /\//, $_;
	    $_="../spcfiles/".$modelname.".".$tmplines[1];
	}

	$lines[$lineIndex]=$_;


	if($_=~/ parameter for the station/){
	    $lineMoments=$lineIndex-1;
	}
	$lineIndex++;
    }

    # multiply the absolute value of moment (in Nm) to the 6 moment "dimensionless" tensors 
    my @tmplines=split ' ', $lines[$lineMoments];
    for($ii=0;$ii<@tmplines;$ii++){
	#print "$ii, $tmplines[$ii]\n";
       

	    
    }
    $tmplines[14]=~s/^e//;
    $tmplines[14]+=7; # Nm -> dyne cm
    $tmplines[14]-=25; # since the DSM units are in 1.e25 dyne cm
    #print "$tmplines[14]\n";
    for($ii=0;$ii<6;$ii++){
	$tmplines[$ii]*=$tmplines[13];
	$tmplines[$ii]=$tmplines[$ii]."e".$tmplines[14];
	#print "$ii, $tmplines[$ii]\n";
    }

    
    $lines[$lineMoments]=$tmplines[0]." ".$tmplines[1]." ".$tmplines[2]." ".$tmplines[3]." ".$tmplines[4]." ".$tmplines[5]." mt (Mrr, Mrt, Mrp, Mtt, Mtp, Mpp) 1.e25 dyne cm";
    
    #print "$lines[$lineMoments]\n";
    #exit;

    

    $lineIndex=0;

    foreach(@lines){
	print FILE1 "$lines[$lineIndex]\n";
		
	$lineIndex++;


    }
    #exit;

    close FILE1;

}



foreach (@shinffiles){
    #print "this is shinffiles : $_\n";
    $modelname=$_; 
    open FILE, $_;
    local @lines =<FILE>;
    close FILE;
    
    $modelname=~s/\s*$//;

    $modelname=~s/^SH_//;

    $modelname=~s/\.inf$//;

    #mkdir ("../spcfiles/$_", 0755) unless (-e "../spcfiles/$_");
    #$tishfile="../inffiles/".$_.".sh.inf";
    #$parfile="../inffiles/".$_.".par.inf";
   
    $tishfile="../inffiles4malbec/$_";
    $_=$tishfile;
    
    open FILE1, ">$tishfile";

       
    my $lineStart=0;
    my $lineEnd=0;
    my $lineIndex=0;
    my $lineMoments=0;
    foreach(@lines){
	if($lineIndex==17){
	    $_="c";
	}
	$_ =~s/\s*$//;
	#$_ =~s/^..\/spcresult\/model\/1/..\/spcfiles/;
	#$_ =~s/SH.spc$/sh.spc/;
	#print "$_\n";
	if($_=~/.spc$/){
	    my @tmplines=split /\//, $_;
	    $_="../spcfiles/".$modelname.".".$tmplines[1];
	}

	$lines[$lineIndex]=$_;


	if($_=~/ parameter for the station/){
	    $lineMoments=$lineIndex-1;
	}
	$lineIndex++;
    }

    # multiply the absolute value of moment (in Nm) to the 6 moment "dimensionless" tensors 
    my @tmplines=split ' ', $lines[$lineMoments];
    for($ii=0;$ii<@tmplines;$ii++){
	#print "$ii, $tmplines[$ii]\n";
       

	    
    }
    $tmplines[14]=~s/^e//;
    $tmplines[14]+=7; # Nm -> dyne cm
    $tmplines[14]-=25; # since the DSM units are in 1.e25 dyne cm
    #print "$tmplines[14]\n";
    for($ii=0;$ii<6;$ii++){
	$tmplines[$ii]*=$tmplines[13];
	$tmplines[$ii]=$tmplines[$ii]."e".$tmplines[14];
	#print "$ii, $tmplines[$ii]\n";
    }

    
    $lines[$lineMoments]=$tmplines[0]." ".$tmplines[1]." ".$tmplines[2]." ".$tmplines[3]." ".$tmplines[4]." ".$tmplines[5]." mt (Mrr, Mrt, Mrp, Mtt, Mtp, Mpp) 1.e25 dyne cm";
    
    #print "$lines[$lineMoments]\n";
    #exit;

    

    $lineIndex=0;

    foreach(@lines){
	print FILE1 "$lines[$lineIndex]\n";
		
	$lineIndex++;


    }
    #exit;

    close FILE1;
   
}



# make conf files for slurm job submissions



my $numberPSV=0;
my $numberSH=0;
my $numberTasks=8; # all the nodes will calculate 8 jobs for PSV and 8 jobs for SH
my $taskIndex=0;

my $numberScripts=@psvinffiles/$numberTasks+1; # number of slurm conf and slurm files that we should write


my $tipsv="/home/franken/dsm4malbec/tipsv4malbec";
my $tish ="/home/franken/dsm4malbec/tish4malbec";


for($taskIndex=0; $taskIndex<$numberScripts; $taskIndex++){

    my $slurmconf="../submissionfiles/".$taskIndex.".conf";
    my $slurmfile="../submissionfiles/".$taskIndex.".slurm";


    my $header=<<EOF;
#!/bin/sh
#SBATCH --job-name selfie
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=16
#SBATCH --partition cpuall,cpunormal,data
#SBATCH --exclusive
module purge
module load intel/compiler intel/mkl slurm
EOF

    open SLURMINF, ">$slurmfile";
    print SLURMINF "$header";
    print SLURMINF "srun -n 16 -l --multi-prog $slurmconf";
    close SLURMINF;
    open SLURMCONF, ">$slurmconf";

    for($ii=0;$ii<$numberTasks;$ii++){
	if($taskIndex*$numberTasks+$ii<@psvinffiles){
	    print SLURMCONF "$ii $tipsv $psvinffiles[$taskIndex*$numberTasks+$ii] ../trash/workPSV$taskIndex".".$ii\n";
	}
    }
    #print "shinffiles $shinffiles[2]\n";
    
    for($ii=0;$ii<$numberTasks;$ii++){
	$jj=$ii+$numberTasks;
	if($taskIndex*$numberTasks+$ii<@shinffiles){
	    print SLURMCONF "$jj $tish $shinffiles[$taskIndex*$numberTasks+$ii] ../trash/workSH$taskIndex".".$ii\n";
	}
    }

    close SLURMCONF;
    

    # !!!! The command below will submit the script immediately!!!
    system("sbatch $slurmfile"); 

}

exit;
    
