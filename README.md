# AFManalysisMatlab
This repository contains a Matlab script to analyze the AFM data to obtain viscoelastic constants.
These scripts were developed for the following publication:
 "Multi-scale study of the architecture, topography, and mechanics of the
human ovary from prepuberty to menopause: a blueprint for next-generation bioengineering and diagnosis."
 Ouni et al., currently under review in Nature communication.
These files allow to analyze AFM data acquired on JPK AFM microscope ('.OUT' file extension) to calculate the viscoelastic response of tissue from AFM measurements.
It fits double exponential plus linear trend to the data and recovers two viscoelastic time constants.
ForceTimeSpectroscopy_cycle_n.m permits to open one '.OUT' file at a time, and analyze multiple force cycles.
ForceTimeSpetroscopy_cycle_1_batch.m   permits to open multiple '.OUT' files at a time, and analyze only one force cycles.
To run these scripts, you need to install the Matlab programming environment first. The untreated raw data can be accessed on the following Mendeley repository: http://dx.doi.org/10.17632/4mcb999nch.1
More detailed instruction to run the data analysis is provided inside the Matlab scripts. The scripts were written in Matlab2019b version (Windows 10), and other versions were not tested.
Figure 1 explains which data points have to be chosen interactively by the user during the data analysis.
