% Finite difference method for solving system of reaction-diffusion
% equations on a disk
% used symbols are from old notation (RD - cytoplasmic Rac1, RT - membrane-bound active Rac1, G - cytoplasmic GAP, Gc - membrane-bound Rac1-GAP complex)

function model_2D_disk
    % Initialize simulation parameters and spatial domain
    P = initializeParameters();

    % Set up the initial conditions
    ic = initializeConditions(P);

    % Solve the ODE system using ode15s (or ode23s)
    [t, u] = ode15s(@vektor_du, P.tspan, ic, [], P);

    
    % Analyze and display the results
    analyzeAndPlotResults(t, u, P);
        
end

function P = initializeParameters() 
    % Model parameters and simulation settings
    
    P.RAD = 6; % radius
    P.n = 6; % Number of radial grid points inside the cytoplasm
    P.m = 40; % Number of angular grid points
    
    P.h = P.RAD/P.n;  % distance between adjacent radial points (radial step)
    P.dphi = 2*pi/P.m; % angular step
    
    P.phi_points = linspace(0, P.m-1, P.m)*P.dphi + P.dphi/2;
    
    for i = 1:P.n
        P.r(i) = (2*i-1)*P.h/2;
        P.r_plus(i) = P.r(i) + P.h/2;
        P.r_minus(i) = P.r(i) - P.h/2;
    end
    
    P.L = P.n*P.h*2*pi;
    
    % N cytoplasmic and M membrane-bound species
    N = 2;
    M = 2;
    P.points = linspace(1, N*P.n*P.m + M*P.m, N*P.n*P.m + M*P.m);
    
    P.tfinal = 4000; % Time domain length (seconds)
    P.tspan = [0, P.tfinal]; % Time domain   
    
    % Periodic boundary conditions (diffusion terms)
    P.xm1 = [P.m,1:P.m-1];
    P.xp1 = [2:P.m,1];
    
    P = defineModelParameters(P);
end


function P = defineModelParameters(P)
    % Model parameters (reaction and diffusion constants)
    P.Rmax = 200;
    
    P.dR = 30;   
    P.dG = 8;
    
    P.k1 = 1.6e-2;  
    P.k2 = 2.5; 
    P.k3 = 0.3; 
    P.k11 = 3.2e-03;   
    P.k12 = 4.8e-03;

    return;
end

function ic = initializeConditions(P)
    % Initializing the spatial domain 
    x = P.n*P.h*P.phi_points;
    
    % Initial conditions
    a = 0.03;
    RD0 = 0.6*ones(1,P.m*P.n);
    G0 = ones(1,P.m*P.n);
    RT0 = 200*(1+a*cos(2*pi*x/P.L));
    Gc0 = 50*(1+a*cos(2*pi*x/P.L+pi/7.7)); 
      
    ic = [RD0,G0,RT0,Gc0];
    
end

function du = vektor_du(t, u, P)
    % ODE system
    RD = u(1:P.n*P.m);
    RD = reshape(RD, P.m, P.n);
    
    G = u(P.n*P.m+1:2*P.n*P.m); 
    G = reshape(G, P.m, P.n);
    
    RT = u(2*P.n*P.m+1:2*P.n*P.m+P.m);    
    Gc = u(2*P.n*P.m+P.m+1:2*P.n*P.m+2*P.m);
        

    MR = 1 - RT/P.Rmax;
    

    % Set of ODEs
    du = zeros(length(P.points),1);
    
    % Cytoplasm
    for i = 1:P.n-1
        vR = linspace(1, P.m, P.m) + (i-1)*P.m;
        vG = linspace(P.n*P.m+1, P.m*(P.n+1), P.m) + (i-1)*P.m;
       
        if i == 1 
            du(vR) = P.dR/P.h^2/P.r(i)*(P.r_plus(i)*(RD(:,i+1) - RD(:,i))) + ...
                P.dR/P.r(i)^2/P.dphi^2*(RD(P.xm1,i) - 2*RD(:,i) + RD(P.xp1,i));
            du(vG) = P.dG/P.h^2/P.r(i)*(P.r_plus(i)*(G(:,i+1) - G(:,i))) + ...
                P.dG/P.r(i)^2/P.dphi^2*(G(P.xm1,i) - 2*G(:,i) + G(P.xp1,i));
        else
            du(vR) = P.dR/P.h^2/P.r(i)*(P.r_plus(i)*(RD(:,i+1) - RD(:,i)) - P.r_minus(i)*(RD(:,i) - RD(:,i-1))) + ...
                P.dR/P.r(i)^2/P.dphi^2*(RD(P.xm1,i) - 2*RD(:,i) + RD(P.xp1,i));
            du(vG) = P.dG/P.h^2/P.r(i)*(P.r_plus(i)*(G(:,i+1) - G(:,i)) - P.r_minus(i)*(G(:,i) - G(:,i-1))) + ...
                P.dG/P.r(i)^2/P.dphi^2*(G(P.xm1,i) - 2*G(:,i) + G(P.xp1,i));
        end
    end
    vR = linspace(1, P.m, P.m) + (P.n-1)*P.m;
    vG = linspace(P.n*P.m+1, P.m*(P.n+1), P.m) + (P.n-1)*P.m;
    du(vR) = -P.dR/P.h^2/P.r(P.n)*P.r_minus(P.n)*(RD(:,P.n) - RD(:,P.n-1)) - ...
            P.r_plus(P.n)/P.r(P.n)/P.h*(RD(:,P.n).*MR.*(P.k1 + P.k11*RT + P.k12*Gc) - P.k3*Gc) + ...
            P.dR/P.r(P.n)^2/P.dphi^2*(RD(P.xm1,P.n) - 2*RD(:,P.n) + RD(P.xp1,P.n));
    du(vG) = -P.dG/P.h^2/P.r(P.n)*P.r_minus(P.n)*(G(:,P.n) - G(:,P.n-1)) - ...
            P.r_plus(P.n)/P.r(P.n)/P.h*(P.k2*RT.*G(:,P.n) - P.k3*Gc) + ...
            P.dG/P.r(P.n)^2/P.dphi^2*(G(P.xm1,P.n) - 2*G(:,P.n) + G(P.xp1,P.n));
    
    % Membrane
    du(2*P.n*P.m+1:2*P.n*P.m+P.m) = - P.k2*RT.*G(:,P.n) + RD(:,P.n).*MR.*(P.k1+P.k11*RT+P.k12*Gc);
    du(2*P.n*P.m+P.m+1:2*P.n*P.m+2*P.m) = P.k2*RT.*G(:,P.n) - P.k3*Gc;
        
        
end

function analyzeAndPlotResults(t, u, P, tb, te)
    % Processing and visualizing the simulation results

    % Extending the spatial domain for periodic boundary conditions in plots
    phi = linspace(-pi, pi, P.m);
    
     
    RD_1 = u(:,1:P.m);
    G_1 = u(:,P.n*P.m+1:P.n*P.m+P.m);
    RT_1 = u(:,2*P.n*P.m+1:2*P.n*P.m+P.m);
    Gc_1 = u(:,2*P.n*P.m+P.m+1:2*P.n*P.m+2*P.m);
    

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
    ylabel(c,'R_D (\mum^-^1)')
    set(c,'FontSize',14)
    
    figure(3);
    surf(phi,t,RT_1,'edgecolor','none')
    xlim([phi(1) phi(end)])
    enhancePlot();
    c = colorbar;
    ylabel(c,'R_D (\mum^-^1)')
    set(c,'FontSize',14)
    
    figure(44);
    surf(phi,t,Gc_1,'edgecolor','none')
    xlim([phi(1) phi(end)])
    enhancePlot();
    c = colorbar;
    ylabel(c,'R_D (\mum^-^1)')
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


