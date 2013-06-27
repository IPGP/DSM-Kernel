#!/usr/bin/perl

require "./messages.pl";


# Initializing procedures
local $ans1;
&messages::messages("welcome");
&messages::question1(\$ans1);
print "$$ans1 b\n";
if($$ans1==1){
    print "Installing DSM Kernel Suite...\n"; 
    
}   





 
