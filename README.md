In the supplementary materials, we include MATLAB codes and accompanying MATLAB data files (Rac1-DGAP1 Fluorescence Data) that enable the reproduction of the results and images presented in this manuscript. Additionally, we supply movies and microscopy files used for generating Figures 1, 2, S1, S2 and S7.

Rac1-DGAP1 Fluorescence Dataset consists of MATLAB data files obtained by analyzing Rac1* and DGAP1# fluorescence data using Quimp software. Each file corresponds to the analysis of an individual cell. The files are organized into three separate folders: 
1. rot_osc_ordered: for cells that exhibited rotating and oscillating patterns,
2. pol_ordered: for cells with stationary (polar) patterns,
3. random: for cells that did not exhibit any particular pattern.

Below is a concise description of the codes provided:

**model_FDM_1D.m**

Solves a set of partial differential equations (PDEs) by discretizing the spatial domain into 100 equally spaced points. The results include kymographs for all seven dependent variables, an autocorrelogram for R_T (generated using the function autocorrelogramFunction.m, which is provided as a separate file), and the phase portrait showing the relationship between membrane-bound Rac1 and DGAP1. Additionally, this code can also be used to generate Figures 4A, 4B, 4C, and the right part of Fig5A.

**model_FDM_2D.m**

Two-dimensional version of the finite difference code used in model_FDM_1D.m. This code solves the set of equations corresponding to cycle I in the model, using a disk as the spatial domain. It can be used to generate Fig9C.

**LSA_FigS4A.m**

Performs linear stability analysis for a specified set of parameters. The analysis begins by computing the system's steady state or fixed point (using the m_SS.m function provided as a separate file). It then proceeds to calculate the Jacobian matrix J and the eigenvalues of the matrix M_q=J-q^2 D, from which the stability of steady state can be determined. Additionally, it plots the growth rate as a function of wave number, corresponding to FigS5A. A slightly modified version can be used to reproduce the left part of Fig5A, and Fig5B.

**kymographs_exp.m**

This script reads Quimp analysis data (files of the form data_FigX.mat). It then interpolates and smooths the data, and plots kymographs and autocorrelograms for the specified time interval. It is used to generate the experimental parts of Fig6, FigS4, and Fig7.

**PCA.m**

Reads Quimp analysis data from the file data_PCA.mat, interpolates this data across a defined spatial and temporal grid, and proceeds with performing principal component analysis. Following the analysis, it plots elements of Fig8, and Fig10.

**stat_ana_FigS5B.m**

Performs the Kruskal-Wallis statistical test to compare the measured circumferences of cells exhibiting first-order oscillatory patterns, second-order oscillatory patterns, and first-order stationary patterns. The data used for this comparison is contained in the file data_stats_mod_v_size.mat.

**noise_analysis_S6Fig.m**

This MATLAB code analyzes the noise profile of YFP intensity data corresponding to the Rac1* signal. It performs a variance-stabilizing transformation, applies a high-pass Butterworth filter to remove low-frequency patterns, and computes the power spectral density (PSD) of the filtered data. The code then visualizes both the averaged PSD (FigS6) and the time-resolved PSD.

**crosscorrelation_Fig2D_Fig2E.m**

This MATLAB code processes fluorescence data from YFP and mRFP channels corresponding to Rac1* and DGAP1# signals, respectively. It interpolates and smooths the data, then calculates the Pearson correlation coefficients for the first 120 seconds of activity across multiple cells. The results for each cell are stored, and a box plot of the Pearson correlation coefficients is generated to visualize the data (Fig2D). A slight modification can be used to generate Fig2E.


