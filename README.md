# ROME2025_CD_Mechanism

# Simulated millennial-scale climate variability driven by a convection-advection oscillator

## Abstract

The last glacial period, between around 115 and 12 thousand years before present, exhibited strong millennial-scale climate variability. This includes abrupt transitions between cold and warm climates, known as Dansgaard-Oeschger (D-O) cycles. D-O cycles have been linked to switches in dynamical regimes of the Atlantic Overturning Meridional Circulation (AMOC), but the exact mechanisms behind abrupt climate changes and AMOC regime shifts remain poorly understood.

This paper introduces the convection-advection oscillator mechanism to explain the millennial-scale oscillations observed in a set of HadCM3 general circulation model simulations forced with snapshots of deglacial meltwater history. The oscillator can be separated into two components acting on different time scales. The fast convection component responds to changes in vertical stratification in the North Atlantic by activating or deactivating deep water formation sites. The slow advection component regulates the accumulation and depletion of salinity in the North Atlantic.

This oscillator mechanism is triggered under specific background conditions and freshwater release patterns. The freshwater perturbation causes an instability that triggers a global salt reorganisation, modifying the North Atlantic stratification. For a given forcing pattern, the system oscillates if the salt transport can lead to an alternating reactivation and deactivation of the AMOC. Otherwise, the climate settles in a warm or cold steady state. This mechanism expands existing theories of millennial-scale variability and provides a general framework for understanding abrupt climate change in general circulation models.


## Reproducibility & Instructions
You will need to download the data from the CEDA archive accessible following this link https://catalogue.ceda.ac.uk/uuid/1f4ebb2944ec43a39ce6c69a8f1942fb/ to run the scripts. The data should be placed in the empty *data* folder conserving the existing hierarchy (i.e. the data folder should contain the folders *inputs*, *tfgbi* etc.). The data files contains the script inputs and the model outputs. If needed, the location of the data can be modified by changing the "data_folder" variable in the notebooks.

You will need to ensure you have the necessary pacakges to run the scripts. Please create a conda/mamba environment loading the rome25_mechanism.yml environment file. You should not have to install any other package. Please make sure that you run the notebooks from the scripts folder so that the references are correct. Please uncomment the lines after each figures to save them in the *Figures* folder. If you want to reproduce the figures exactly as they are in the manuscript, you will need to have the Inter font (https://fonts.google.com/specimen/Inter) installed on your machine.

If you experience an issue with this repository, please log an issue on the corresponding GitHub page (https://github.com/Olnavy/ROME2024_CD_Mechanism/).

## Experiment names

| Experiment name | Experiment labels | Comments |
| ----------- | ----------- | ----------- | 
| XOUPA | CTRL | Initial 1000 years corresponding to spin-up |
| XOUPH | 21k | - |
| TFGBI | 20.7k | - |
| XOUPF | 18.2k | - |
| XOUPK | 20.7_tdc | - |
| XPPBF | 20.7_no_gst | - |

## References
No references to declare.
