function [xp, vp, rho_back, Q, QM] = initilization_Two_Stream_Basic(params, N)
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
% Perturbation in particle positions
xp=xp + A*(L/N)*sin(2*pi*xp/L*mode);

% Maxwellian velocity 
vp=VT*randn(N,1);
pm=[1:N]';
pm=1-2*mod(pm,2);
% add the beam drift velocity
vp=vp+pm.*V0;

end