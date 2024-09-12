%% Define the directory containing the data files
myFolder = 'Insert folder path'; % Define your working folder, e.g. C:\My Drive\Image analysis\Variables\
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
end
filePattern = fullfile(myFolder, '*.mat');
matFiles = dir(filePattern);


results = struct(); % Initialize results structure
allPearsonCoefficients = []; % Initialize an array for all Pearson coefficients


%% Process each file
for k = 1:length(matFiles)
    baseFileName = matFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    matData = load(fullFileName);  % Load each data file

    % Extract time interval information and prepare time vector
    delta_t = matData.ans.FI;
    if strcmp(baseFileName, '2014_4_29_Series011.mat')
        delta_t = 4;
    end
    tq = delta_t * linspace(0, length(matData.ans.fluo)-2, length(matData.ans.fluo)-1);
    
    
    % Median contour length from Quimp analysis
    perim = round(median(matData.ans.stats(:,8)));
    if strcmp(baseFileName, 'R&D_181127_s006_35-406.mat')
        perim = 44;
    end

    %% Interpolation 

    numInterpPoints = 100;  % Number of interpolation points along x
    xx = linspace(-2 * perim, 3 * perim, 5 * numInterpPoints);

    yfp_E_int = zeros(length(tq), length(xx));
    mrfp_E_int = zeros(size(yfp_E_int));

    for i = 1:length(tq)
        x = perim * matData.ans.outlines{i,1}(:,1)';
        k0 = find(x == 0);
        x = circshift(x, -k0 + 1);  % Aligning spatial coordinates (start at x = 0)

        MeanCytoFluo_YFP = matData.ans.fluoStats(i,1,7); % Mean cytoplasmic fluorescence
        MeanCytoFluo_RFP = matData.ans.fluoStats(i,2,7);
        yfp = circshift(matData.ans.fluo{i,1}(:,1,1)', -k0 + 1)*MeanCytoFluo_YFP; % use *MeanCytoFluo_YFP if data is normalized to interior 
        mrfp = circshift(matData.ans.fluo{i,1}(:,2,1)', -k0 + 1)*MeanCytoFluo_RFP;
        
        if strcmp(baseFileName, 'R&D_190521_006_1.mat') || strcmp(baseFileName, 'R&D_190527_001_2.mat') || strcmp(baseFileName, 'R&D_190528_006_2.mat') || strcmp(baseFileName, 'R&D_190528_008_1.mat') % in these cases there was a slight leakage of yfp signal (Rac1) into mrfp channel, so subtracted roughly 3% of yfp signal from mrfp signal (based on the emission spectrum of yfp) 
            c = 0.03;

            % Calculate the initial vector mrfp1
            mrfp1 = mrfp - c * yfp;

            % Loop to check and adjust p until there are no negative values in mrfp1
            while any(mrfp1 < 0)
                % Lower the value of p by 0.005
                c = c - 0.005;

                % Recalculate mrfp1 with the new value of p
                mrfp1 = mrfp - c * yfp;
            end
            cc(i) = c;
            mrfp = mrfp1;
        end

        % Extending data points for more accurate smoothing
        xl = x - perim; xr = x + perim; 
        xl2 = x - 2 * perim; xr2 = x + 2 * perim;
        x_P = [xl2, xl, x, xr, xr2];
        yfp_P = repmat(yfp, 1, 5); 
        mrfp_P = repmat(mrfp, 1, 5); 

        % Cubic smoothing spline interpolation
        yfp_E_int(i,:) = csaps(x_P, yfp_P, 1, xx); 
        mrfp_E_int(i,:) = csaps(x_P, mrfp_P, 1, xx);
    end
     
  %% Cubic smoothing spline smoothing

    smoothingParam = 0.1;
    xt = {tq, xx}; 

    yfp_P_smooth = csaps(xt, yfp_E_int, smoothingParam, xt);
    mrfp_P_smooth = csaps(xt, mrfp_E_int, smoothingParam, xt);

    % Triming to the original interval
    index0 = find(abs(xx) == min(abs(xx)));
    index1 = find(abs(xx - perim) == min(abs(xx - perim)));
    
    y_int = yfp_E_int(:, index0:index1);
    m_int = mrfp_E_int(:, index0:index1);  

    y_smooth = yfp_P_smooth(:, index0:index1);
    m_smooth = mrfp_P_smooth(:, index0:index1);  
    
    xq = linspace(0, perim, index1 - index0 + 1);
    
    %% Define time segments for Pearson correlation calculation
    dt = 120; % time of each segment (seconds)
    t_n = round(tq(end) / dt);
    t_a = tq(end) / t_n;
    tb_vek = linspace(0, tq(end) - t_a, t_n); % Start times of segments
    te_vek = linspace(0, tq(end) - t_a, t_n) + t_a; % End times of segments

    pearsonCoefficients = zeros(1, length(tb_vek)); % Array to store Pearson coefficients for each segment

    % Calculate Pearson correlation for each segment
    for j = 1:length(tb_vek)
        tb = tb_vek(j);
        te = te_vek(j);
        in_tb = find(tq >= tb, 1, 'first');
        in_te = find(tq <= te, 1, 'last');

        y_a = y_smooth(in_tb:in_te, :);
        m_a = m_smooth(in_tb:in_te, :);

        if ~isempty(y_a) && ~isempty(m_a)
            pearsonCoefficients(j) = corr2(y_a, m_a);
        else
            pearsonCoefficients(j) = NaN; % Handle cases where data might be missing
        end
    end

    % Store results
    results(k).fileName = baseFileName;
    results(k).pearsonCoefficients = pearsonCoefficients;
    
    % Concatenate the current file's Pearson coefficients to the global array
    allPearsonCoefficients = [allPearsonCoefficients, pearsonCoefficients]; % Append current file's coefficients
    
    % First Pearson coefficient (first 120 seconds)
    P_120s(k) = pearsonCoefficients(1);
    

end

med_all_PC = median(allPearsonCoefficients);
perc25_all_PC = prctile(allPearsonCoefficients,25);
perc75_all_PC = prctile(allPearsonCoefficients,75);

%% Plotting
jitterStrength = 0.05; % Adjust this value based on your preference for separation
jitter = jitterStrength * (rand(size(P_120s)) - 0.5); % Uniform jitter around the fixed X value
figure;
scatter(ones(size(P_120s)) + jitter, P_120s, 60, 'filled','k','MarkerEdgeColor', 'w', 'LineWidth', 1); % Apply jitter
ylabel('Mean Pearson coefficient');
xlim([0.9 1.1]); % Adjust the limits appropriately to account for jitter
ylim([-1 1]);
set(gca, 'XTickLabel',[]);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
hold on
% Box plot
h = boxplot(P_120s, 'Positions', 1, 'Widths', 0.07, 'Colors', 'b', 'Symbol', ''); % Overlay the box plot without outliers

% Customize box plot
set(h, 'LineWidth', 1.5); % Set the line width

% Customizing individual elements of the box plot
boxHandles = findobj(gca, 'Tag', 'Box');
medianHandles = findobj(gca, 'Tag', 'Median');
whiskerHandles = findobj(gca, 'Tag', 'Whisker');
capsHandles = findobj(gca, 'Tag', 'Cap');

% Set properties
set(boxHandles, 'Color', 'k', 'LineWidth', 1.5); % Set box color and line width
set(medianHandles, 'Color', 'r', 'LineWidth', 1.5); % Set median line color and width
set(whiskerHandles, 'Color', 'k', 'LineWidth', 1.5); % Set whisker color and line width
set(capsHandles, 'Color', 'k', 'LineWidth', 1.5); % Set caps color and line width
xlim([0.93 1.07])
