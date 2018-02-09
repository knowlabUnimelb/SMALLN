# SMALLN
Simulation Code for Smith &amp; Little - Small is Beautiful: In Defense of the Small-N Design

simulateSmallN.m - code for simulating additive factors model with N = 4
- Running this file will cretae smallNsimulation.mat

plotSmallN.m - load results of smallNsimulation.mat and plot

simulateLargerN.m - code for simulating additive factors model with N > 4
- To match paper, set N = [8, 16, 32, 64, 128]
- Running this file will create largerNsimulation_[N].mat 

plotLargeNsimulations.m - load results of largerNsimulations_[N].mat and plot

-----

Several .mat files have been uploaded. These store the output of the simulation .m files

-----

Additional files needed to run the code:

aggregate.m, allcomb.m, displayTime.m, getbandwidth.m, keep.m, mstrfind.m
