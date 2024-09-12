% Calculating steady state
function F = m_SS(SS,k1,k2,k3,k4,k5,k6a,k6b,k11,k12,k41,k42,Rmax,Dmax,Ru,Gu,Du)

SS1 = SS(1); % RD
SS2 = SS(2); % G
SS3 = SS(3); % D
SS4 = SS(4); % RT
SS5 = SS(5); % Dm
SS6 = SS(6); % Gc
SS7 = SS(7); % Dc

% Steady state needs to satisfy the following equations
F(1) = k3*SS6 + k6b*SS7 - SS1*(1-SS4/Rmax)*(k1 + k11*SS4 + k12*SS6);
F(2) = k3*SS6 - k2*SS2*SS4;
F(3) = -(1-SS5/Dmax)*SS3*(k4+k41*SS4+k42*SS6) + (k6a+k6b)*SS7;
F(4) = k5*SS4*SS5 - (k6a+k6b)*SS7;
F(5) = SS1 + SS4 + SS6 + SS7 - Ru;
F(6) = SS2 + SS6 - Gu;
F(7) = SS3 + SS5 + SS7 - Du;
end

