function damprate = calc_Landau_Maxwell(time, histEnergy, dx)

	histQOI = Energy_Norms.fourierSpectrum(histEnergy, 2);
    % Compute critical points and damping rate of energies
    lnE      = log(histQOI);
    lnE_big  = [lnE(1), lnE, lnE(end)];
    peaks    = find(diff(lnE_big(1:end-1)) > 0 & diff(lnE_big(2:end)) < 0);
    Epeaks   = lnE(peaks);
    Tpeaks   = time(peaks);
    damprate = -(Epeaks(4) - Epeaks(3))/(Tpeaks(4) - Tpeaks(3));
end

