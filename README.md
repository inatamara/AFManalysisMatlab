# AFManalysisMatlab
This repository contains Matlab script to analyse the AFM data to obtain elastic, and viscoelastic constants.
These scripts were developped for the following publication:
 "Multi-scale study of the architecture, topography and mechanics of the
human ovary from prepuberty to menopause: a blueprint for next-generation bioengineering and diagnosis"
 Ouni et al., currently under review in Nature communication.
The file ForceTimeSpetroscopy_AFMviscoelasticModel.m allows to analyse AFM data aquired on JPK AFM microscope ('.OUT' file extension) to
calculate the viscoelastic response of a tissue from AFM measurements.
It fits double exponential plus linear trend to the data and recovers two viscoelastic time constants.
To run this scripts you need to install Matlab programming enviroment first. The analysed data and raw untreated data can be accessed on the following Mendeley repository: http://dx.doi.org/10.17632/4mcb999nch.1
More detailed instruction to run the data analysis is provided inside the matlab scripts.
The scripts were written in Matlab2019b version (Windows 10), and other versions were not tested.
The PlotResultsViscoelatic3ageGroups.m permits to plot the final results that can be accessed in the following Mendeley data repository: http://dx.doi.org/10.17632/4mcb999nch.1
