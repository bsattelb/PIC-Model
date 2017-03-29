clear;
clc;

addpath('..')

% params = [L, A, VT, V0]
params = [4*pi, 0.5, 0, 0.2];
paramNames = {'L', 'A', 'VT', 'V0'};
QOI = 'Growth Rate';
change = [true, true, false, true];
Nparams = sum(change);
alpha = 0.05;
max_vals = (1+alpha)*params;
min_vals = (1-alpha)*params;

DT = 0.01;
NT = 20/DT+1;
NG = 32;
N = 1e3;

PIC = PIC.PIC_setup(DT, NT, NG, N, 'Two_Stream_Basic');

h = 1e-6;
Nsamples = 100;

[w, output, Xs, sdev, Atrials] = Sensitivity.Linear_Fit(max_vals, min_vals, Nsamples, PIC, 'test_params', change);

Sensitivity.plotter_Linear_Fit(Nparams, Nsamples, paramNames, QOI, w, output, Xs)

save('Results_LF/Two_Stream_Basic.mat')

rmpath('..')