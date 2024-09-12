%% Define the directory containing the data files
myFolder = 'Insert folder path'; % Define your working folder, e.g. C:\My Drive\Image analysis\Variables\
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
end
filePattern = fullfile(myFolder, '*.mat');
matFiles = dir(filePattern);


PSD_y_int_filtered_tot = [];
%% Process each file
for k = 1  % calculates noise profile for one cell - if calculating statistics for batch of cells then use k = 1:length(matFiles)
    baseFileName = matFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    matData = load(fullFileName);  % Load each data file

    % Extract time interval information and prepare time vector
    delta_t = matData.ans.FI;
    if strcmp(baseFileName, 'R&D_140429_011_0_350.mat')
        delta_t = 4;
    end
    tq = delta_t * linspace(0, length(matData.ans.fluo)-2, length(matData.ans.fluo)-1);
    
    
    % Median contour length from Quimp analysis
    perim = round(median(matData.ans.stats(:,8)));
    if strcmp(baseFileName, 'R&D_181127_s006_35-406.mat')
        perim = 44;
    end

    %% Interpolation using csaps

    numInterpPoints = 100;  % Number of interpolation points along x
    xx = linspace(-2 * perim, 3 * perim, 5 * numInterpPoints);

    yfp_E_int = zeros(length(tq), length(xx));
    mrfp_E_int = zeros(size(yfp_E_int));

    for i = 1:length(tq)
        x = perim * matData.ans.outlines{i,1}(:,1)';
        k0 = find(x == 0);
        x = circshift(x, -k0 + 1);  % Aligning spatial coordinates (start at x = 0)

        MeanCytoFluo_YFP = matData.ans.fluoStats(i,1,7); % Mean cytoplasmic fluorescence
        yfp = circshift(matData.ans.fluo{i,1}(:,1,1)', -k0 + 1)*MeanCytoFluo_YFP; % use *MeanCytoFluo_YFP if data is normalized to interior 
        
        % Extending data points for more accurate smoothing
        xl = x - perim; xr = x + perim; 
        xl2 = x - 2 * perim; xr2 = x + 2 * perim;
        x_P = [xl2, xl, x, xr, xr2];
        yfp_P = repmat(yfp, 1, 5); 

        % Cubic smoothing spline interpolation
        yfp_E_int(i,:) = csaps(x_P, yfp_P, 1, xx); 
    end
    
    % Triming to the original interval
    index0 = find(abs(xx) == min(abs(xx)));
    index1 = find(abs(xx - perim) == min(abs(xx - perim)));
    
    xq = linspace(0, perim, index1 - index0 + 1);
    
    y_int = yfp_E_int(:, index0:index1);
    m_int = mrfp_E_int(:, index0:index1); 
     

%% noise analysis 
% Parameters
N_time = size(y_int, 1);  % Number of time points
N_space = size(y_int, 2); % Number of spatial points
Fs_space = 100 / perim;      % Spatial sampling frequency (samples per unit length). If done for one cell at a time, then Fs = 100 / perim. If analysis id done for the batch of cells, then Fs_space = 2 is a good approximation

% Step 1: Variance-Stabilizing Transformation (Square Root Transformation)
y_int_sqrt = sqrt(y_int);  % Apply square root transformation to stabilize variance

% Step 2: Design a High-Pass Filter to Remove Low-Frequency Patterns
Fc = 0.15;  % Cutoff frequency in cycles per unit length (adjust based on your data)
[b, a] = butter(4, Fc/(Fs_space/2), 'high');  % 4th order Butterworth high-pass filter

% Step 3: Apply the High-Pass Filter to the Transformed Data
y_int_filtered = zeros(size(y_int_sqrt));  % Preallocate for filtered data

for t = 1:N_time
    y_int_filtered(t, :) = filtfilt(b, a, y_int_sqrt(t, :));  % Apply zero-phase filter to each time point
end

% Step 4: Compute the Power Spectral Density (PSD) of the Filtered Data Over Time
freq_space = 0:(Fs_space/N_space):(Fs_space/2);  % Frequency axis
PSD_y_int_filtered = zeros(N_time, N_space/2+1);  % Preallocate for PSD

for t = 1:N_time
    Y_filtered = fft(y_int_filtered(t, :));  % Compute FFT for each time point
    PSD_y_int_filtered(t, :) = 2*abs(Y_filtered(1:N_space/2+1)).^2 / N_space;  % Compute single-sided PSD
end

PSD_y_int_filtered_tot = [PSD_y_int_filtered_tot;PSD_y_int_filtered]; % use for batch of cells 

% Step 5: Average the PSD Across Time to Focus on Overall Noise Characteristics
avg_PSD_filtered = mean(PSD_y_int_filtered, 1);

% Plot the Averaged Power Spectral Density of the Filtered Data
figure;
plot(freq_space, avg_PSD_filtered,'LineWidth', 2, 'Color', 'k');
xlabel('Spatial Frequency (cycles per micrometer)');
ylabel('Power/Frequency');
grid off;

% Time-Resolved PSD (Heatmap)
% Visualize how the PSD changes over time using a heatmap
figure;
imagesc(1:N_time, freq_space, PSD_y_int_filtered');
axis xy;  % Flip the y-axis for correct frequency orientation
colorbar;
title('Time-Resolved Power Spectral Density of Noise');
xlabel('Time');
ylabel('Spatial Frequency (cycles per unit length)');
set(gca, 'YScale', 'linear');  % Optional: You can set to 'log' for better frequency resolution

si = size(y_int_filtered);
% Visualize the filtered data (after removing low-frequency patterns)
figure;
surf(xq,linspace(0, 1, si(1)),real(y_int_filtered),'edgecolor','none')
colorbar;
title('Filtered Fluorescence Intensity Data (Noise)');
xlabel('Space');
ylabel('Time');
view(2)



end

%% Use this for batch of cells
% avg_PSD_filtered_tot = mean(PSD_y_int_filtered_tot, 1);
% 
% figure (7);
% plot(freq_space, avg_PSD_filtered_tot,'LineWidth', 2, 'Color', 'k');
% xlabel('Spatial Frequency (cycles per micrometer)');
% ylabel('Power/Frequency');
% grid off;
