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
% Lorentzian velocity distribution - inverse cdf sampling
%vp=VT*tan(pi*(rand(N, 1) - 0.5)) + V0;
%vp=VT*trnd(1,N,1) + V0;
vp=VT*(randn(N,1)./randn(N,1))+V0;

%remove outliers
temp = abs(vp) > 5;
while nnz(temp) > 0
	vp(temp) = VT*(randn(nnz(temp),1)./randn(nnz(temp),1))+V0;
	temp = abs(vp) > 5;
end

% Perturbation in particle positions
xp=xp + A*cos(2*pi*xp/L*mode);

end