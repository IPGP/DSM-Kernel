#!/bin/bash
mpirun -np 8 /home/dtp02/fuji/DSM-KS2011//bin/mpiSGTpsv < /home/dtp02/fuji/DSM-KS2011//inf_calcul//input10.inf > /home/dtp02/fuji/DSM-KS2011//log/10.ak135.psv.log
mpirun -np 8 /home/dtp02/fuji/DSM-KS2011//bin/mpiSGTsh  < /home/dtp02/fuji/DSM-KS2011//inf_calcul//input10.inf > /home/dtp02/fuji/DSM-KS2011//log/10.ak135.sh.log
