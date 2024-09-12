%% Kymographs and autocorrelograms 
clear all

%% Load Quimp analysis data
varName = 'Insert file name';  % e.g. data_Fig4A
path1 = 'Insert folder path'; 
data = load(strcat(path1,varName));  

% Time resolution
delta_t = data.ans.FI;

% Median contour length from Quimp analysis
perim = round(median(data.ans.stats(:,8)));

%% Interpolation 

% Time and space domains for interpolation
tq = delta_t * linspace(0, length(data.ans.fluo) - 2, length(data.ans.fluo) - 1);

numInterpPoints = 100;  % Number of interpolation points along x
xx = linspace(-2 * perim, 3 * perim, 5 * numInterpPoints);

yfp_E_int = zeros(length(tq), length(xx));
mrfp_E_int = zeros(size(yfp_E_int));

for i = 1:length(tq)
    x = perim * data.ans.outlines{i,1}(:,1)';
    k0 = find(x == 0);
    x = circshift(x, -k0 + 1);  % Aligning spatial coordinates (start at x = 0)
    
    MeanCytoFluo_YFP = data.ans.fluoStats(i,1,7); % Mean cytoplasmic fluorescence
    MeanCytoFluo_RFP = data.ans.fluoStats(i,2,7);
    yfp = circshift(data.ans.fluo{i,1}(:,1,1)', -k0 + 1)*MeanCytoFluo_YFP; % use *MeanCytoFluo_YFP if data is normalized to interior 
    mrfp = circshift(data.ans.fluo{i,1}(:,2,1)', -k0 + 1)*MeanCytoFluo_RFP;
    
    % Extending data points for more accurate smoothing
    xl = x - perim; xr = x + perim; 
    xl2 = x - 2 * perim; xr2 = x + 2 * perim;
    x_P = [xl2, xl, x, xr, xr2];
    yfp_P = repmat(yfp, 1, 5); 
    mrfp_P = repmat(mrfp, 1, 5); 
    
    % Cubic smoothing spline interpolation
    yfp_E_int(i,:) = csaps(x_P, yfp_P, 1, xx); 
    inn1 = find(yfp_E_int>110);
    mrfp_E_int(i,:) = csaps(x_P, mrfp_P, 1, xx);
end

%% Cubic smoothing spline smoothing

smoothingParam = 0.5;
xt = {tq, xx}; 

yfp_P_smooth = csaps(xt, yfp_E_int, smoothingParam, xt);
mrfp_P_smooth = csaps(xt, mrfp_E_int, smoothingParam, xt);

% Triming to the original interval
index0 = find(abs(xx) == min(abs(xx)));
index1 = find(abs(xx - perim) == min(abs(xx - perim)));

y_smooth = yfp_P_smooth(:, index0:index1);
m_smooth = mrfp_P_smooth(:, index0:index1);

xq = linspace(0, perim, index1 - index0 + 1);

%%

font = 14;

% Plot y_smooth
figure(1);
surf(xq, tq, y_smooth, 'edgecolor', 'none')
set(gca, 'FontSize', font)
xlim([xq(1) xq(end)])
ylim([0 tq(end)])
xlabel('x(\mum)')
ylabel('t(s)')
title('DPAKa(GBD)-DYFP')
view(2)
c = colorbar;
ylabel(c, 'Intensity (AU)')

% Plot m_smooth if there is mRFP channel (test based on its uniformity)
if m_smooth(1, 1) ~= mean(mean(m_smooth))
    figure(2);
    surf(xq, tq, m_smooth, 'edgecolor', 'none')
    set(gca, 'FontSize', 14)
    xlim([xq(1) xq(end)])
    ylim([0 tq(end)])
    xlabel('x(\mum)')
    ylabel('t(s)')
    title('DGAP1-mRFP')
    view(2)
    c = colorbar;
    ylabel(c, 'Intensity (AU)')
end



%% Autocorrelation calculation
tb = 150; te = 470;  % Start and end times for autocorrelation calculation
in_tb = find(abs(tq - tb) == min(abs(tq - tb)), 1);
in_te = find(abs(tq - te) == min(abs(tq - te)), 1);
y_smooth_a = y_smooth(in_tb:in_te, :);

[m, n] = size(y_smooth_a); 

y_mean_x = mean(y_smooth_a, 2);  % Average across each row (x variable)
dy = y_smooth_a - y_mean_x;  % Deviations from the mean for each time point

% Squared deviations and mean square deviation
dy_2 = dy.^2; 
avg_dy_2 = mean(dy_2, 'all'); 

% Autocorrelation calculation
Cy = xcorr2(dy) / avg_dy_2 / numel(y_smooth_a);

dxq = xq(2)-xq(1);
dtq = tq(2)-tq(1);

% New axes for plotting the autocorrelation
x1 = linspace(-pi,pi,n);  % šetanje po pola x-osi
t1 = linspace(-tq(m)/2,tq(m)/2,m);
Cy_half = Cy(1+(m-1)/2:m+(m-1)/2,1+(n-1)/2:n+(n-1)/2);

figure(3);
phi = linspace(-pi,pi,100);
surf(phi,tq(in_tb:in_te)-tq(in_tb),y_smooth_a,'edgecolor','none')
set(gca,'FontSize',font)
xlim([phi(1) phi(end)])
ylim([0 tq(in_te)-tq(in_tb)])
title('DPAKa(GBD)-DYFP (subinterval)')
xlabel('\phi (rad)')
ylabel('t(s)')
c = colorbar;
ylabel(c, 'Intensity (AU)')
set(c,'FontSize',font)
view(2)

figure(4);
surf(x1,t1,Cy_half,'edgecolor','none')
set(gca,'FontSize',font)
view(2)
xlabel('\Delta\phi (rad)')
ylabel('\Deltat (s)')
xlim([x1(1) x1(end)]); ylim([t1(1) t1(end)]);
c = colorbar;
set(c,'FontSize',font)
ylabel(c, 'Autocorrelation')
colormap jet
set(c,'FontSize',font)


