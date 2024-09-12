function autocorrelogramFunction(u, t, n, tb, te)

    % Interpolating to a regular grid (since time-steps in original data aren't uniform)
    tq = linspace(0, t(end), numel(t));
    
    % Membrane-bound Rac1
    RT = horzcat(u(:, 3*n+1:4*n), u(:, 3*n+1)); 
    siz = size(RT);
    
    RT_interp = zeros(numel(tq), siz(2));
    for j = 1:siz(2)
        RT_interp(:, j) = interp1(t, RT(:, j), tq, 'linear');
    end

    in_tb = find(abs(tq - tb) == min(abs(tq - tb)), 1);
    in_te = find(abs(tq - te) == min(abs(tq - te)), 1);

    RT_a = RT_interp(in_tb:in_te, :);

    % Normalization
    [m, n] = size(RT_a);
    y_mean_x = mean(RT_a, 2); % Vectorized mean calculation
    dy = RT_a - y_mean_x; % Difference from the mean
    dy_2 = dy.^2;
    avg_dy_2 = mean(dy_2, 'all'); 

    Cy = xcorr2(dy) / avg_dy_2 / n / m; % Autocorrelation

    x1 = linspace(-pi, pi, n);
    t1 = linspace(-tq(m)/2, tq(m)/2, m);
    phi = linspace(-pi, pi, length(x1));
    t_plot_a = tq(in_tb:in_te) - tq(in_tb);

    % Visualization
    for figNum = 1:2
        figure(figNum)

        if figNum == 1
            surfData = RT_a;
            xData = phi;
            yData = t_plot_a;
        else
            surfData = Cy(1+(m-1)/2:m+(m-1)/2, 1+(n-1)/2:n+(n-1)/2);
            xData = phi;
            yData = t1;

        end

        surf(xData, yData, surfData, 'edgecolor', 'none');
        set(gca, 'FontSize', 14, 'TickDir', 'out');
        view(2);
        xlabel('\phi (rad)');
        colorbar('FontSize', 14);
        xlim([phi(1) phi(end)]);
        ylim([yData(1), yData(end)]);

        if figNum == 1
            ylabel('t(s)')
            c = colorbar;
            ylabel(c,'R_T (\mum^-^1)')
        else
            ylabel('t(s)')
            colormap jet;
            ylabel('\Deltat (s)')
            c = colorbar;
            ylabel(c,'Autocorrelation')
        end
    end
end
