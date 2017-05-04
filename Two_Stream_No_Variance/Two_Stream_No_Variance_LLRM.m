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
NT = 70/DT+1;
NG = 128;
N = 1e4;

PIC = PIC.PIC_setup(DT, NT, NG, N, 'Two_Stream_Basic');

Nsamples = 10000;
Nsamples2 = 2500;
p = 5000;

[evalues, U, output, Xs, Xs2, graddamp, sdev, Atrials] = ...
    Sensitivity.Local_Linear(max_vals, min_vals, Nsamples, ...
                                   Nsamples2, p, PIC, 'test_params', change);

Sensitivity.plotter_Local_Linear(Nparams, Nsamples, paramNames, QOI, evalues, U, output, Xs)

save('Results_LLRM/Two_Stream_Basic.mat')

rmpath('..')