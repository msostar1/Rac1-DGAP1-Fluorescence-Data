function model_FDM
    % Initialize simulation parameters and spatial domain
    P = initializeParameters();

    % Set up the initial conditions
    ic = initializeConditions(P);

    % Solve the ODE system using ode15s (or ode23s)
    [t, u] = ode15s(@vektor_du, P.tspan, ic, [], P);

    % Define the time subinterval for autocorrelation analysis and phase portrait 
    % (te must be samller than P.tfinal)
    tb = 1500;  % Beginning of the subinterval
    te = 1850;  % End of the subinterval
    
    % Check if te is larger than P.tfinal and warn if needed
    if te > P.tfinal
        warning('The value of te (%.2f) is larger than P.tfinal (%.2f). Please adjust it to be smaller.', te, P.tfinal);
    else
        % Analyze and display the results
        analyzeAndPlotResults(t, u, P, tb, te);
        % Call the autocorrelogram function if te is within a valid range
        autocorrelogramFunction(u, t, P.nx, tb, te);
    end
end

function P = initializeParameters() 
    % Model parameters and simulation settings
    P.tfinal = 2000; % Time domain length (seconds)
    P.L = 40; % Spatial domain length
    P.nx = 100; % Number of grid points
    
    P.tspan = [0, P.tfinal]; % Time domain
    P.deltax = P.L / P.nx; % Space step
    
    % Indices for the spatial domain 
    P.x = 1:P.nx;
    P.x1 = P.x;
    P.x2 = P.x + P.nx;  
    P.x3 = P.x + 2*P.nx;
    P.x4 = P.x + 3*P.nx;
    P.x5 = P.x + 4*P.nx;  
    P.x6 = P.x + 5*P.nx;
    P.x7 = P.x + 6*P.nx;
    
    % Periodic boundary conditions (diffusion terms)
    P.xm1 = [P.nx,1:P.nx-1];
    P.xp1 = [2:P.nx,1];
    
    P = defineModelParameters(P);
end


function P = defineModelParameters(P)
    % Model parameters (reaction and diffusion constants)
    P.Rmax = 200;
    P.Dmax = 40;

    P.dR = 30;   
    P.dG = 8;
    P.dD = 12;

    P.k1 = 2e-2;  
    P.k2 = 0.22; 
    P.k3 = 0.45; 
    P.k11 = 4e-03;   
    P.k12 = 6e-03;

    P.k5 = 0.095;
    P.k6a = 0.6;
    P.k6b = 0.43;
    P.k4 = 0.24;

    P.k41 = 0.8e-03;
    P.k42 = 3.6e-03;
    return;
end

function ic = initializeConditions(P)
    % Initializing the spatial domain 
    x = linspace(-P.L / 2, P.L / 2 - P.L / P.nx, P.nx);
    
    % Initial conditions
    m = 1;
    a = 0.03;
    RD0 = 130*(1+a*cos(m*2*pi*x/P.L)); RT0 = 0*(1+a*cos(m*2*pi*x/P.L));
    Gc0 = 0*(1+a*cos(m*2*pi*x/P.L));  G0 = 30*(1+a*cos(m*2*pi*x/P.L));
    Dc0 = 0*(1+a*cos(m*2*pi*x/P.L)); D0 = 40*(1+a*cos(m*2*pi*x/P.L));
    Dm0 = 0*(1+a*cos(m*2*pi*x/P.L));

    ic = [RD0,G0,D0,RT0,Dm0,Gc0,Dc0];
end

function du = vektor_du(t, u, P)
    % ODE system
    RD = u(P.x1);
    G = u(P.x2);  
    D = u(P.x3);
    RT = u(P.x4);
    Dm = u(P.x5);
    Gc = u(P.x6);
    Dc = u(P.x7);

    MR = 1 - RT/P.Rmax;
    MD = 1 - Dm/P.Dmax;

    % Set of ODEs
    du = zeros(7*P.nx,1);

    du(P.x1) = P.k6b*Dc + P.k3*Gc - RD.*MR.*(P.k1+P.k11*RT+P.k12*Gc) + P.dR/P.deltax^2*(RD(P.xm1)-2*RD(P.x)+RD(P.xp1));
    du(P.x2) = P.k3*Gc - P.k2*RT.*G + P.dG/P.deltax^2*(G(P.xm1)-2*G(P.x)+G(P.xp1));
    du(P.x3) = (P.k6a + P.k6b)*Dc - MD.*D.*(P.k4+P.k41*RT + P.k42*Gc) + P.dD/P.deltax^2*(D(P.xm1)-2*D(P.x)+D(P.xp1));
    du(P.x4) = P.k6a*Dc - P.k2*RT.*G - P.k5*RT.*Dm + RD.*MR.*(P.k1+P.k11*RT+P.k12*Gc);
    du(P.x5) = MD.*D.*(P.k4+P.k41*RT + P.k42*Gc) - P.k5*RT.*Dm;
    du(P.x6) = P.k2*RT.*G - P.k3*Gc;
    du(P.x7) = P.k5*RT.*Dm - (P.k6a + P.k6b)*Dc;
end

function analyzeAndPlotResults(t, u, P, tb, te)
    % Processing and visualizing the simulation results

    % Extending the spatial domain for periodic boundary conditions in plots
    phi = linspace(-pi, pi, P.nx+1);
    
    RD_1 = horzcat(u(:,1:P.nx),u(:,1));
    G_1 = horzcat(u(:,P.nx+1:2*P.nx),u(:,P.nx+1));  
    D_1 = horzcat(u(:,2*P.nx+1:3*P.nx),u(:,2*P.nx+1));
    RT_1 = horzcat(u(:,3*P.nx+1:4*P.nx),u(:,3*P.nx+1));
    Dm_1 = horzcat(u(:,4*P.nx+1:5*P.nx),u(:,4*P.nx+1));
    Gc_1 = horzcat(u(:,5*P.nx+1:6*P.nx),u(:,5*P.nx+1));
    Dc_1 = horzcat(u(:,6*P.nx+1:7*P.nx),u(:,6*P.nx+1));

    % Plots
    figure(1);
    surf(phi,t,RD_1,'edgecolor','none')
    xlim([phi(1) phi(end)])
    enhancePlot();
    c = colorbar;
    ylabel(c,'R_D (\mum^-^1)')
    set(c,'FontSize',14)
    
    figure(2);
    surf(phi,t,G_1,'edgecolor','none')
    xlim([phi(1) phi(end)])
    enhancePlot();
    c = colorbar;
    ylabel(c,'G (\mum^-^1)')
    set(c,'FontSize',14)
    
    figure(3);
    surf(phi,t,D_1,'edgecolor','none')
    xlim([phi(1) phi(end)])
    enhancePlot();
    c = colorbar;
    ylabel(c,'D (\mum^-^1)')
    set(c,'FontSize',14)
    
    figure(4);
    surf(phi,t,RT_1,'edgecolor','none')
    xlim([phi(1) phi(end)])
    enhancePlot();
    c = colorbar;
    ylabel(c,'R_T (\mum^-^1)')
    set(c,'FontSize',14)
    
    figure(5);
    surf(phi,t,Dm_1,'edgecolor','none')
    xlim([phi(1) phi(end)])
    enhancePlot();
    c = colorbar;
    ylabel(c,'D_m (\mum^-^1)')
    set(c,'FontSize',14)
    
    figure(6);
    surf(phi,t,Gc_1,'edgecolor','none')
    xlim([phi(1) phi(end)])
    enhancePlot();
    c = colorbar;
    ylabel(c,'G_c (\mum^-^1)')
    set(c,'FontSize',14)
    
    figure(7);
    surf(phi,t,Dc_1,'edgecolor','none')
    xlim([phi(1) phi(end)])
    enhancePlot();
    c = colorbar;
    ylabel(c,'D_c (\mum^-^1)')
    set(c,'FontSize',14)
    
    % Plotting the phase portrait
    in_tb = find(abs(t-tb)==min(abs(t-tb))); 
    in_te = find(abs(t-te)==min(abs(t-te)));
    
    RT_phase = RT_1(in_tb:in_te,1)';
    D_phase = Dm_1(in_tb:in_te,1)' + Dc_1(in_tb:in_te,1)';
    z = ones(size(RT_phase));
    col = linspace(0,te-tb,length(RT_phase));  % Color - varies with time
    figure(8)
    colormap jet
    surface([RT_phase;RT_phase],[D_phase;D_phase],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2,'linesmoothing', 'on');
        
    x_lim_1 = min(RT_phase) - 0.1*(max(RT_phase)-min(RT_phase)); 
    x_lim_2 = max(RT_phase) + 0.1*(max(RT_phase)-min(RT_phase));
    y_lim_1 = min(D_phase) - 0.1*(max(D_phase)-min(D_phase)); 
    y_lim_2 = max(D_phase) + 0.1*(max(D_phase)-min(D_phase));
    xlim([x_lim_1 x_lim_2]); 
    ylim([y_lim_1 y_lim_2]);
    xlabel('R_T (\mum^-^1)')
    ylabel('D_m + D_c (\mum^-^1)')
    set(gca,'FontSize',14)
    set(gca,'TickDir','in')
    c = colorbar;
    ylabel(c,'t (s)')
    set(c,'FontSize',14)
end

% Plot parameters
function enhancePlot()
    set(gca,'FontSize',14)
    view(2)
    set(gca,'TickDir','out')
    xticks([-2 0 2])
    xticklabels({'-2','0','2'})
    xlabel('\phi (rad)')
    ylabel('t(s)')
    grid off;
end


