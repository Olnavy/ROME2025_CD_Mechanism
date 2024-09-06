# ROME2024_CD_Mechanism

# Simulated millennial-scale climate variability driven by a convection-advection oscillator

## Abstract

The last glacial period, between $\sim$ 115 and 12 thousand years before present, exhibited strong millennial-scale climate variability. This includes abrupt transitions between cold and warm climates, known as Dansgaard-Oeschger (D-O) events. D-O events have been linked to switches in regimes of the Atlantic Overturning Meridional Circulation (AMOC), but the exact mechanisms behind abrupt climate changes and AMOC regimes switches remain poorly understood.

This paper introduces the convection-advection oscillator mechanism to explain the millennial-scale oscillations observed in a set of HadCM3 general circulation model simulations forced with snapshots of deglacial meltwater history. The oscillator can be separated into two components acting on different time scales. The fast convection component responds to changes in vertical stratification in the North Atlantic by activating or deactivating its deep water formation sites. The slow advection component regulates the accumulation and depletion of salinity in the North Atlantic.

This oscillator mechanism is triggered under specific background conditions and freshwater release combinations. The initial perturbation introduces an instability that triggers a global salt reorganisation, modifying the North Atlantic stratification. For a given configuration, the system oscillates if the salt redistribution can lead to both the reactivation and the deactivation of the AMOC. Otherwise, the climate settles in a warm or cold steady state. This new mechanism expands the existing millennial-scale variability theories and provides a general framework for understanding abrupt climate changes in general circulation models.


## Reproducibility & Instructions
Please refer to the paper to download the data necessary to run the notebooks. The data are accessible on the CEDA archive.

The Figures can be reproduced by running the notebooks after loading the rome2024_mechanism environment in conda.

If you want to reproduce exactly the Figures, you will need to have the Inter font (https://fonts.google.com/specimen/Inter) installed.

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
No references
