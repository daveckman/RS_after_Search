"Reusing search data in ranking and selection: What could possibly go wrong?"


Documentation
-------------
The included MATLAB scripts are used to test simulation experiments and generate figures for paper. The three main files are makePlots.m, plotRealSearchg, and RealSearch.m


makePlots.m
-----------
Tests Bechhofer, Rinott, Modified Gupta, and Screen-to-the-Best on Adversarial Search (AS) and plots Figures 1 and 2.
Calls subroutines such as AdvSearch and AdvSearch Rinott that return data from AS for k systems.


plotRealSearchg.m
-----------------
Plots Figure 3(a).


RealSearch.m
------------
Tests Modified Gupta on a realistic search and plots Figure 3(b).
Calls RealSearchLog as a subroutine that returns data from realistic search for k systems.

Other files
-----------
calcBechhoferh.m: Calculate h_B (Bechhofer h) as defined in paper
Bech_h_list.txt: Text file with list of calculated Bechhofer h for settings of k used in experiments.
calcRinott.m: Calculate h_R (Rinott h) as defined in paper
calcDblInt.m: Calculate double integral in expression for Rinott h. Subroutine of calcRinott.m
Rinott_h_list.txt: Text file with list of calculated Rinott h for settings of k used in experiments.