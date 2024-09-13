Oscillatory dynamics of Rac1 activity in Dictyostelium discoideum amoebae

Authors: Marko Šoštar, Maja Marinović, Vedrana Filić, Nenad Pavin and Igor Weber

**Supporting information**

Folder “Movies”

S1 Movie. Representative rotating monopole. A Rac1*-enriched domain travels along the cell membrane (corresponds to Fig 1A).

S2 Movie. Representative rotating dipole. Two distinct Rac1*-enriched domains travel along the cell membrane (corresponds to S1 Fig A).

S3 Movie. Representative oscillating monopole. A Rac1*-enriched domain periodically relocates from one side of the cell to the opposite side (corresponds to Fig 1B).

S4 Movie. Representative oscillating dipole. Two distinct Rac1*-enriched domains periodically undergo 90-degree orientation shifts (corresponds to S1 Fig B).

S5 Movie. Representative stationary monopole. A Rac1*-enriched domain remains on one side of a migrating cell (corresponds to Fig 1C).

S6 Movie. Representative stationary dipole. Two distinct Rac1*-enriched domains located on opposite sides of the cell maintain their positions (corresponds to S1 Fig C).

S7 Movie. Representative rotating monopoles. Rac1* and DGAP1-enriched domains, located on opposite sides of the cell, travel along the cell membrane (corresponds to Fig 2A). 

S8 Movie. Representative oscillating dipoles. Two distinct Rac1*-enriched domains and oppositely localized DGAP1-enriched domains periodically shift their orientation by 90 degrees (corresponds to S2 Fig).

S9 Movie. Representative stationary monopoles out of phase. Rac1*-enriched domain and oppositely localized DGAP1-enriched domain maintain their positions (corresponds to Fig 2B).

S10 Movie. Representative stationary monopoles in phase. Rac1*-enriched and DGAP1-enriched domains are co-localized, maintaining their position on the cell membrane (corresponds to Fig 2C).



Folder “Supplementary Figures”

S1 Fig. Representative bipolar patterns of cortical domains enriched in Rac1*. 

S2 Fig. A representative oscillating dipole in a double-labeled cell.

S3 Fig. Impact of diffusion coefficients on pattern formation in the Rac1-GAP system.

S4 Fig. Comparison of the experimentally observed patterns of Rac1* with the patterns of Rac1T derived from computer simulations.

S5 Fig. Effect of circumference size on pattern formation.

S6 Figure. Representative power spectral density of the Rac1 activity.

S7 Fig. Monopolar patterns of Rac1T in DGAP1/GAPA knock-out cells. 


Folder “Matlab Codes” 

Contains all MATLAB codes used to generate data shown in Figures 3, 4, 5, 6, S3 and S4.

In the following, a concise description of the MATLAB codes is provided:

model_FDM_1D.m

Solves a set of partial differential equations (PDEs) by discretizing the spatial domain into 100 equally spaced points. The results include kymographs for all seven dependent variables, an autocorrelogram for R_T (generated using the function autocorrelogramFunction.m, which is provided as a separate file), and the phase portrait showing the relationship between membrane-bound Rac1 and DGAP1. Additionally, this code can also be used to generate Figures 4A, 4B, 4C, and the right part of Fig5A.

model_FDM_2D.m

Two-dimensional version of the finite difference code used in model_FDM_1D.m. This code solves the set of equations corresponding to cycle I in the model, using a disk as the spatial domain. It can be used to generate Fig9C.

LSA_FigS4A.m

Performs linear stability analysis for a specified set of parameters. The analysis begins by computing the system's steady state or fixed point (using the m_SS.m function provided as a separate file). It then proceeds to calculate the Jacobian matrix J and the eigenvalues of the matrix M_q=J-q^2 D, from which the stability of steady state can be determined. Additionally, it plots the growth rate as a function of wave number, corresponding to FigS5A. A slightly modified version can be used to reproduce the left part of Fig5A, and Fig5B.

kymographs_exp.m

This script reads Quimp analysis data (files of the form data_FigX.mat). It then interpolates and smooths the data, and plots kymographs and autocorrelograms for the specified time interval. It is used to generate the experimental parts of Fig6, FigS4, and Fig7.

PCA.m

Reads Quimp analysis data from the file data_PCA.mat, interpolates this data across a defined spatial and temporal grid, and proceeds with performing principal component analysis. Following the analysis, it plots elements of Fig8, and Fig10.

stat_ana_FigS5B.m

Performs the Kruskal-Wallis statistical test to compare the measured circumferences of cells exhibiting first-order oscillatory patterns, second-order oscillatory patterns, and first-order stationary patterns. The data used for this comparison is contained in the file data_stats_mod_v_size.mat.

noise_analysis_S6Fig.m

This MATLAB code analyzes the noise profile of YFP intensity data corresponding to the Rac1* signal. It performs a variance-stabilizing transformation, applies a high-pass Butterworth filter to remove low-frequency patterns, and computes the power spectral density (PSD) of the filtered data. The code then visualizes both the averaged PSD (FigS6) and the time-resolved PSD.
crosscorrelation_Fig2D_Fig2E.m

This MATLAB code processes fluorescence data from YFP and mRFP channels corresponding to Rac1* and DGAP1# signals, respectively. It interpolates and smooths the data, then calculates the Pearson correlation coefficients for the first 120 seconds of activity across multiple cells. The results for each cell are stored, and a box plot of the Pearson correlation coefficients is generated to visualize the data (Fig2D). A slight modification can be used to generate Fig2E.


Folder “Matlab Data”: Data obtained from the QUIMP analysis of the cell outline fluorescence intensity used as the input for MATLAB routines to generate Figures 4, 5, 6, S3 and S4.

Folder “Microscopy Data”: Confocal microscopy series used to generate Figures 1, S1, 2, S2, and S7

Folder “Raw Microscopy Data”: Raw microscopy data used to generate Figures 1, S1, 2, S2, and S7

Folder “pol_ordered”: Matlab data files obtained from the QUIMP analysis of cells exhibiting stationary polarization, used for generating Fig. 2D.

Folder “rot_osc_ordered”: Matlab data files obtained from the QUIMP analysis of cells exhibiting rotating and oscillating patterns, used for generating Fig. 2D.

Folder “random”: Matlab data files obtained from the QUIMP analysis of cells exhibiting random behavior.
