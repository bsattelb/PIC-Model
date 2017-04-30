clear;
clc;

addpath('..')

% params = [L, A, VT, V0]
params = [4*pi, 0.5, 1, 0];
paramNames = {'L', 'A', 'VT', 'V0'};
QOI = 'Damping Rate';
change = [true, true, true, true];
Nparams = sum(change);
alpha = 0.05;
max_vals = (1+alpha)*params;
min_vals = (1-alpha)*params;
max_vals(end) = 0.01;
min_vals(end) = -0.01;

DT = 0.025;
NT = 10/DT;
NG = 128;
N = 1e6;

PIC = PIC.PIC_setup(DT, NT, NG, N, 'Landau_Maxwell');

Nsamples = 20;
Nsamples2 = 15;
p = 17;

[evalues, U, output, Xs, Xs2, graddamp, sdev, Atrials] = ...
    Sensitivity.Local_Linear_Model(max_vals, min_vals, Nsamples, ...
                                   Nsamples2, p, PIC, 'Averaging', true, ...
                                   'Atolerance', 5e-4);


Sensitivity.plotter_Active_Subspaces(Nparams, Nsamples, paramNames, QOI, evalues, U, output, Xs)

save('Results_AS/Landau_Maxwell.mat')

rmpath('..')