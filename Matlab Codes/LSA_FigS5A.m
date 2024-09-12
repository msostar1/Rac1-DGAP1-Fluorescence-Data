% Linear stability analysis (fig. S4A)

%% Model parameters

Rmax = 200;
Dmax = 40;

dR = 30;   
dG = 8;
dD = 12;

k1 = 2e-2;  
k2 = 0.22; 
k3 = 0.45; 
k11 = 4e-03;   
k12 = 6e-03;

k5 = 0.095;
k6a = 0.6;
k6b = 0.43;
k4 = 0.24;

k41 = 0.8e-03;
k42 = 3.6e-03;

R_tot = 130; % Total Rac1 linear concentration  
G_tot = 30; % Total GAP linear concentration  
D_tot = 40; % Total DGAP1 linear concentration  

%% Steady state

ini_g = repmat(50, 7, 1); % Initial guess
equationSystem = @(SS) m_SS(SS, k1, k2, k3, k4, k5, k6a, k6b, k11, k12, k41, k42, Rmax, Dmax, R_tot, G_tot, D_tot);

% Finding a non-negative solution
while true
    SS = fsolve(equationSystem, ini_g);
    if min(SS) >= 0
        break;  % Exit the loop if all elements are non-negative
    end
    ini_g = 100 * rand(7, 1);  % Generate a new starting guess
end

%% Jacobian matrix

syms U1 U2 U3 U4 U5 U6 U7 k1s k2s k3s k4s k5s k6as k6bs k11s k12s k41s k42s Rmaxs Dmaxs
% jacobian(f,U) computes the Jacobian matrix of symbolic function f with respect to U
JS = jacobian([k3s*U6 + k6bs*U7 - U1*(1-U4/Rmaxs)*(k1s + k11s*U4 + k12s*U6),...    
    k3s*U6 - k2s*U2*U4,...
    -U3*(1-U5/Dmaxs)*(k4s + k41s*U4 + k42s*U6) + (k6as + k6bs)*U7,...
    k6as*U7 - k2s*U2*U4 - k5s*U4*U5 + U1*(1-U4/Rmaxs)*(k1s + k11s*U4 + k12s*U6),...
    U3*(1-U5/Dmaxs)*(k4s + k41s*U4 + k42s*U6) - k5s*U4*U5,...
    -k3s*U6 + k2s*U2*U4,...
    k5s*U4*U5 - (k6as + k6bs)*U7], [U1, U2, U3, U4, U5, U6, U7]);

% Converting a symbolic expression to a function
J = matlabFunction(JS, 'Vars', [U1,U2,U3,U4,U5,U6,U7,k1s,k2s,k3s,k4s,k5s,k6as,k6bs,k11s,k12s,k41s,k42s,Rmaxs,Dmaxs]);

% Evaluating the Jacobian matrix (in steady state)
J_SS = J(SS(1),SS(2),SS(3),SS(4),SS(5),SS(6),SS(7),k1, k2, k3, k4, k5, k6a, k6b, k11, k12, k41, k42, Rmax, Dmax);

%% Stability of steady state

L = logspace(log10(pi), 3, 2000); % Spatial domain length (cell cirumference) 
q = 2 * pi ./ L; % Corresponding wave numbers

s = size(J_SS,1); 
n = numel(q);
M = zeros(s,s,n);
lambda = zeros(s, n);
lam1 = zeros(n, 1);

% D as a 7x7 diffusion matrix with dR, dG, dD on its diagonal
D = diag([dR, dG, dD, zeros(1, s-3)]);

% Iteratively adjusting the Jacobian matrix and calculating eigenvalues
for i = 1:n   
    M(:,:,i) = J_SS - q(i)^2 * D;
    lambda(:, i) = eig(M(:,:,i));  % Eigenvalues
    
    % Maximum real parts of eigenvalues
    lam1(i) = max(real(lambda(:, i)));
end

%% Plot

figure (1)
plot(q,lam1)
set(gca,'FontSize',14) 
xlabel('q (\mum^-^1)')
ylabel('\sigma_q (s^-^1)')

