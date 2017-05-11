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

Nsamples = 15;

[w, output, Xs, sdev, Atrials] = Sensitivity.Linear_Fit( ...
    max_vals, min_vals, Nsamples, PIC, 'test_params', change, ...
    'Averaging', true, 'Atolerance', 5e-3);

Sensitivity.plotter_Linear_Fit(Nparams, Nsamples, paramNames, QOI, w, output, Xs, sdev)

save('Results_LF/Landau_Maxwell.mat')

rmpath('..')