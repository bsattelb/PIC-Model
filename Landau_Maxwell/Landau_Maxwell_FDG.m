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

h = 1e-2;
Nsamples = 15;

[evalues, U, output, outputplus, Xs, graddamp, sdev, Atrials] = ...
    Sensitivity.FD_Gradient(max_vals, min_vals, h, Nsamples, PIC, ...
    'test_params', change, 'Averaging', true, 'Atolerance', 1e-3);

Sensitivity.plotter_FD_Gradient(Nparams, Nsamples, paramNames, QOI, evalues, U, output, Xs)

save('Results_AS/Landau_Maxwell.mat')

rmpath('..')