package messages;

sub messages{
    my ($type) = @_;
    
    
    if($type=="welcome"){
	print ("##########################################################################\n");
	print ("#                                                                        #\n");
	    
	    print ("#       DSM Kernel Suite 2012 (beta)                                     #\n");
	    print ("#                                                                        #\n");
	    print ("#                                                                        #\n");
	    print ("#                                        2012.2. FUJI Nobuaki            #\n");
	    print ("#                                                                        #\n");
	    print ("#                           after programs by Kenji Kawai & Li Zhao      #\n");
	    print ("#                                                                        #\n");
	    print ("#                           with incredible helps by                     #\n");
	    print ("#                           Robert J. Geller, Dimitri Komatitsch,        #\n");
	    print ("#                           Sebastien Chevrot, Vadim Monteiller,         #\n");
	    print ("#                           Hiromitsu Mizutani, Nozomu Takeuchi,         #\n");
	    print ("#                           Kensuke Konishi, Roland Martin & Marie Calvet#\n");
	    print ("#                                                                        #\n");
	    print ("##########################################################################\n");
	print ("#\n");
	print ("#\n");
	print ("# DSM Kernel Suite consists of:\n");
	print (" #0 Synthetic Calculation (simple version): tish,tipsv,spcsac\n");
	print (" #1 Green Function Calculation: mpiSGTpsv/sh\n");
	    print (" #2 Kernel Calculation: mpiKernelMaker\n");
	    print (" #3 Visualisation: ReadKernel/plotker\n");
	    print ("#\n");
	print ("# DSM Kernel Suite requires:\n");
	print ("#   mandatory: (mpi)f90, perl\n");
	print ("#   desired: Matlab, GMT\n");
	    
	print ("##########################################################################\n");
	}
}

sub question1{
    (*stadin) =@_;

    print ("\n1. Installation\n2. Parameter Configuration\n3. Calculation\n");
    print ("Please indicate what you'd like to do:\n");
    $$stadin =<STDIN>;
    $$stadin =~s/\s*$//;
    print "$$stadin a\n";

    
}




1;
