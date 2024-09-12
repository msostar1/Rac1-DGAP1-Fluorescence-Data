% Principal component analysis
clear all

%% Define variables and paths
varName = 'data_PCA'; % File name
path1 = 'Insert folder path'; 
dataFile = fullfile(path1, [varName, '.mat']); % Use fullfile for path construction

% Load data
loadedData = load(dataFile);  
delta_t = loadedData.ans.FI; % Frame interval (seconds)

%% Interpolation on an extended x-axis
perim = round(median(loadedData.ans.stats(:,8))); % Median cell contour length
t_out = delta_t * linspace(0, numel(loadedData.ans.fluo) - 2, numel(loadedData.ans.fluo) - 1);
nx = 100; % Number of points for x-axis interpolation
xx = linspace(-2 * perim, 3 * perim, 5 * nx); 

% Initialize matrices for interpolated data
numFrames = numel(loadedData.ans.fluo) - 1;
yfp_P_int = zeros(numFrames, numel(xx));
mrfp_P_int = zeros(numFrames, numel(xx));

% Iterate over each frame for processing
for i = 1:numFrames
    % Expand and replicate spatial coordinates and intensities
    x = perim * circshift(loadedData.ans.outlines{i}(:,1)', [0, -1]);
    yfp = circshift(loadedData.ans.fluo{i}(:,1,1)' * loadedData.ans.fluoStats(i,1,7), [0, -1]);
    x_P = repmat(x, 1, 5) + repelem(perim * (-2:2), numel(x));
    yfp_P = repmat(yfp, 1, 5);

    % Cubic smoothing spline interpolation
    yfp_P_int(i, :) = csaps(x_P, yfp_P, 0.7, xx);

    % Check and process mRFP data if varying
    if loadedData.ans.fluo{1}(1,2,1) ~= mean(loadedData.ans.fluo{1}(:,2,1))
        mrfp = circshift(loadedData.ans.fluo{i}(:,2,1)' * loadedData.ans.fluoStats(i,2,7), [0, -1]);
        mrfp_P = repmat(mrfp, 1, 5);
        mrfp_P_int(i, :) = csaps(x_P, mrfp_P, 0.7, xx);
    end
end

%% Trimming to the original x-interval and selecting the t-interval

tb = 0; te = 820;  % Start and end times
[~, in_tb] = min(abs(t_out - tb)); 
[~, in_te] = min(abs(t_out - te));
[~, index0] = min(abs(xx - 0));  % Find the index of the value in xx closest to 0
[~, index1] = min(abs(xx - perim));  % Find the index of the value in xx closest to perim
yfp_int = yfp_P_int(in_tb:in_te, index0:index1);

% Check for fluorescence data variations and adjust if necessary
if loadedData.ans.fluo{1}(1,2,1) ~= mean(loadedData.ans.fluo{1}(:,2,1))
    mrfp_int = mrfp_P_int(in_tb:in_te, index0:index1);
    % Apply correction for potential leakage into the mrfp channel
    mrfp_int = mrfp_int - 0.04 * yfp_int; 
end

%% Kymographs

t_out = t_out(in_tb:in_te) - t_out(in_tb);
s = size(yfp_int,2);
xq1 = linspace(-pi,pi,s);

font0 = 14;

figure(1);
surf(xq1,t_out,yfp_int,'edgecolor','none');
set(gca,'FontSize',font0)
xlim([xq1(1) xq1(end)])
ylim([0 t_out(end)])
set(gca,'TickDir','out')
xticks([-2 0 2])
xticklabels({'-2','0','2'})
xlabel('\phi (rad)')
ylabel('t(s)')
title('DPAKa(GBD)-DYFP')
c = colorbar;
ylabel(c, 'Intensity (AU)')
set(c,'FontSize',font0)
view(2)

if loadedData.ans.fluo{1}(1,2,1) ~= mean(loadedData.ans.fluo{1}(:,2,1))
    figure(2);
    surf(xq1,t_out,mrfp_int,'edgecolor','none');
    set(gca,'FontSize',font0)
    xlim([xq1(1) xq1(end)])
    ylim([0 t_out(end)])
    set(gca,'TickDir','out')
    xticks([-2 0 2])
    xticklabels({'-2','0','2'})
    xlabel('\phi (rad)')
    ylabel('t(s)')
    title('DGAP1-mRFP')
    c = colorbar;
    ylabel(c, 'Intensity (AU)')
    set(c,'FontSize',font0)
    view(2)
end

%% PCA variances
% Check if there is a significant variation in the mRFP channel (second channel)
if loadedData.ans.fluo{1,1}(1,2,1) ~= mean(loadedData.ans.fluo{1}(:,2,1))
    ym = [yfp_int; mrfp_int];  % Concatenate yfp and mrfp matrices vertically
else
    ym = yfp_int;
end

% Centering and normalizing the matrix ym
ym = (zscore(ym'))';

% Principal Component Analysis with time points as variables
[COEFF,SCORE,latent] = pca(ym');  

% Calculating the percentage of variance explained by each principal component
PC_variances_perc = 100 * latent / sum(latent);  

% Cumulative percentage of explained variance
kum_PC_var_perc = cumsum(PC_variances_perc(1:10));

% Plotting the results
font1 = 14;
figure(3)  
bar(PC_variances_perc(1:10), 0.7, 'FaceColor', [.8 .8 .8], 'EdgeColor', [0 0 0])
set(gca, 'FontSize', font1, 'XLim', [0.4 10.6])
xlabel('Principal Component', 'FontSize', font1)
ylabel('Percentage of Variance Explained', 'FontSize', font1)
title('PCA Analysis: Variance Explained by Each PC', 'FontSize', font1)

%% Determining oscillation periods from the first PC coefficients  
v1_yfp = COEFF(1:length(t_out), 1);  % First PC coefficients for YFP (Rac1)
% Smoothing
p = 0.015;  % Smoothing parameter
v1_yfp_smooth = csaps(t_out, v1_yfp', p, t_out);  

% Autocorrelation for YFP PC1 coefficients and determining periods
[cor_v1_y, lag_v1_y] = xcorr(v1_yfp - median(v1_yfp), 'coeff');
[pks_v1_y, locs_v1_y] = findpeaks(cor_v1_y);
periods_v1_y = delta_t * diff(locs_v1_y);

font = 14;
    
figure (4) % 1st principal component for yfp
plot(t_out,v1_yfp_smooth,'linewidth',1)
set(gca,'FontSize',font)
xlabel('Time (s)')
ylabel('PC no.1 coefficients')
title('Rac1 PC no.1')
xlim([0,820])

figure(5)
plot(delta_t*lag_v1_y,cor_v1_y)
set(gca,'FontSize',font)
xlabel('\Deltat (s)')
ylabel('Autocorrelation')

period_estimate = 260; % estimated from periods_v1_y

%% Averaging Dynamics

% Projection vectors for the first and second principal components
c1 = ym' * COEFF(:, 1);  % Projection of data onto the first PC
c2 = ym' * COEFF(:, 2);  % Projection of data onto the second PC

% Determine the rotation angle to align the phase of the first point with zero
theta = -atan2(c2(1), c1(1));
rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];  % Rotation matrix
C = [c1'; c2'];  % Combine projections into a matrix
CC = rot * C;  % Apply rotation

% Calculate phase relationship among points
phase = atan2(CC(2, :)', CC(1, :)');
phase = phase - phase(1);  % Normalize phase to start from zero

% Phase difference mapped to time index differences
omega = 2 * pi / period_estimate;
d_t = round(-phase / omega / delta_t);

% Extending the time axis for signal translation
extended_time = [t_out(1) - flip(t_out(2:end)), t_out, t_out(end) + (t_out(2:end) - t_out(1))];

% Initialize matrix for translated signal values
N = NaN(nx, length(extended_time));
brt = length(t_out);

N(1,brt+1:2*brt) = yfp_int(:,1); % Central portion of matrix N is set to be equal to the first column of matrix yfp_int 
for j = 2:nx 
    N(j,brt+1-d_t(j):2*brt-d_t(j)) = yfp_int(:,j);
end

% Calculate and clean the average dynamics
N_mean = nanmean(N);
inN = find(isnan(N_mean)); 
N_mean(inN) = []; 
yfp_mean = N_mean;

% Final extended time vector for averaged dynamics
tpr = delta_t * (0:(length(yfp_mean) - 1));

% Smoothing averaged dynamics
py_mean = 0.001;  % Smoothing parameter for yfp_mean
yfp_mean_smooth = csaps(tpr, yfp_mean, py_mean, tpr);

% Handling mRFP channel if there is one
if loadedData.ans.fluo{1}(1,2,1) ~= mean(loadedData.ans.fluo{1}(:,2,1))
    M = NaN([nx,length(extended_time)]);
    M(1,brt+1:2*brt) = mrfp_int(:,1);
    for j = 2:nx 
        M(j,brt+1-d_t(j):2*brt-d_t(j)) = mrfp_int(:,j);
    end
    
    M_mean = nanmean(M);
    inM = find(isnan(M_mean)); 
    M_mean(inM) = []; 
    mrfp_mean = M_mean;
    
    ko = 0.004; % Bleaching correction coefficient (needs to be determined separately)
    pm_mean = 0.00005;  % Smoothing parameter for mrfp_mean
    mrfp_mean_smooth = csaps(tpr, mrfp_mean, pm_mean, tpr) + ko * tpr;
end

figure(6) % Projection of data onto the first two principal components
scatter(CC(1,:),CC(2,:),25,'.','k')
hold on
xline(0, 'LineWidth', 0.5);  % X-axis
yline(0, 'LineWidth', 0.5);  % Y-axis
box on
hold off

% Find columns in N that don't consist entirely of NaN values
notAllNaNs = ~all(isnan(N), 1);
% Find the index of the first and last column where not all entries are NaN
firstColIndex = find(notAllNaNs, 1, 'first');

figure (7) % The aligned time series of the fluorescence intensities for DPAKa(GBD)-DYFP along with averaged time series
plot(extended_time - extended_time(firstColIndex),N,'o','LineStyle','none','markersize',1)
hold on
plot(tpr,yfp_mean_smooth,'linewidth',2,'color',[.2 .2 .2])
hold off

if loadedData.ans.fluo{1}(1,2,1) ~= mean(loadedData.ans.fluo{1}(:,2,1))
    figure(8) % Smoothed average yfp i mrfp time series 
    plot(tpr,yfp_mean_smooth,tpr,mrfp_mean_smooth)
    set(gca,'FontSize',font)
    xlabel('Time (s)')
    ylabel('Intensity (AU)')
    legend('DPAKa(GBD)-DYFP','DGAP1-mRFP')
    
    
    add = repmat(ko*extended_time,numel(nx),1); % mrfp signal bleaching correction
    figure (9) % The aligned time series of the fluorescence intensities for DGAP1-mRFP along with averaged time series
    plot(extended_time - extended_time(firstColIndex),M+add,'o','LineStyle','none','markersize',1)
    hold on
    plot(tpr,mrfp_mean_smooth,'linewidth',2,'color',[.2 .2 .2])
    hold off
    
    
    % Subinterval for phase diagram
    tb = 50; % Begining time
    te = 1060; % End time
    in_tb = find(abs(tpr-tb)==min(abs(tpr-tb))); 
    in_te = find(abs(tpr-te)==min(abs(tpr-te)));
    yfp_phase = yfp_mean_smooth(in_tb:in_te);
    mrfp_phase = mrfp_mean_smooth(in_tb:in_te);
    z2 = zeros(size(yfp_phase));
    x_lim_l = min(yfp_phase) - 0.1*(max(yfp_phase)-min(yfp_phase)); x_lim_r = max(yfp_phase) + 0.1*(max(yfp_phase)-min(1));
    y_lim_d = min(mrfp_phase) - 0.1*(max(mrfp_phase)-min(mrfp_phase)); y_lim_u = max(mrfp_phase) + 0.1*(max(mrfp_phase)-min(mrfp_phase));
    col = linspace(0,tpr(in_te)-tpr(in_tb),length(yfp_phase));  % Color - varies with time
    figure(10) % Phase diagram
    colormap jet
    surface([yfp_phase;yfp_phase],[mrfp_phase;mrfp_phase],[z2;z2],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',1,'linesmoothing', 'on');
        xlim([x_lim_l x_lim_r]); ylim([y_lim_d y_lim_u]);
    set(gca,'FontSize',font)
    box on
    xlabel('DPAKa(GBD)-DYFP')
    ylabel('DGAP1-mRFP')
    c = colorbar;
    ylabel(c,'t (s)') 
end
