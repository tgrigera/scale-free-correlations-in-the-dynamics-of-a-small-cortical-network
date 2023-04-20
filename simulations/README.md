# Code for the Supplementary Information of "Scale-free correlations..."

This directory contains the code used to produce figures 9 to 13 of "Scale-free correlations in the dynamics of a small ($N \sim 10000$) cortical network", by 
S. Camargo et al.

The programs given here are slight modifications of the code in [this Github repository](https://github.com/DanielAlejandroMartin/Kappa_C), which was orinially developed for the work reported in

  - Trejo, E. J. A., Martin, D. A., De Zoysa, D., Bowen, Z., Grigera, T. S., Cannas, S. A., Losert, W., Chialvo, D. R., Finite-size correlation behavior near a critical point: A simple metric for monitoring the state of a neural network. _Phys. Rev. E_ __106,__ 0543143.  (2022). [DOI 10.1103/PhysRevE.106.054313](http://dx.doi.org/10.1103/PhysRevE.106.054313)


Each computation requires several simple codes.

Running the following scripts in the following order will allow to reproduce the results in the above mentioned figures.


1. In "1-Dynamics_Generator", the dynamics for different computations are generated, using different scripts.  "Source" folder contains the program "Program_Dyn.f90". 

   - ScriptMultiRunOriginal: 
	 It is used for Normal, Convolved, Noisy_Convolved (Fig. 9) and Shuffle (Fig. 13). It runs 
	 copies of "../Source/Progran_Dyn.f90"

	- ScriptMultiRun_TB::
	 It is used for Time Binning (Fig. 10). Essentially the same codea as for [Original",
	 but records the summed activity over 20 steps.

	- ScriptMultiRunArousal:: 
	 It is used for Arousal (Fig. 11). Essentially the same codea as for "Original", but changes r_1 every 5000 steps

	- ScriptShuffler::
	 To be run after "ScriptMultiRunOriginal". It generates a shuffled signal from "Original"

2. In 2A-Compute_CCF/CCF1-CorrelationAnalisys:

	- Run ScriptMultiRun<folder> to compute the CCF of that folder where folder may be: "Original" "Shuffle" "Arousal" "TimeBinning" "Convol" or "ConvolNoise" "Convol" and "ConvolNoise" use "Original" simulations as input


3. In 2B-Compute_CCF/CCF2-ZeroFinder:
The average of C(r) over different windows and the zeros crossings of C(r) are computed for each 
floder in CCF1-CorrelationAnalisys
	
	- Run ScriptCCF2A-ScriptJoinWindows, ScriptCCF2B-FindZeros ND ScriptCCF2C-AverageFunctPlots 

4. In 2B-TimeCorrel/AC1-ExtractActivity/

	- Run ScriptMultiRunOriginal and ScriptMultiRunConvol

5. In 2B-TimeCorrel/AC2-Autocorr$ 

	- Run Script2Autocorr

6. In 3-Figures: Figures are ploted with gnuplot.

	- Run gnuplot "GnuFig<number>"


