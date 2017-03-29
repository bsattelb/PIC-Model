function [xp, vp, rho_back, Q, QM] = initilization_Landau_Maxwell(params, N)
% perturbation amplitude and mode
L= params(1);

% simulation box length
A=params(2);

% thermal velocity
VT=params(3);

% beam velocity
V0 = params(4);


% electron charge to mass ratio
QM=1;
% perturbation mode
mode=1;
% computational particle charge
Q=1/(QM*N/L);
%background charge given by background (not moving) ions
rho_back=-Q*N/L;

% uniform particle distribution in space
xp=linspace(0,L-L/N,N)';
% Maxwellian velocity 
vp=VT*randn(N,1) + V0;

% Perturbation in particle positions
xp=xp + A*cos(2*pi*xp/L*mode);

end