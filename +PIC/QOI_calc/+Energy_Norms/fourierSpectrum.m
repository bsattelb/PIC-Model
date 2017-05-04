function histMode = fourierSpectrum(histEnergy, mode)
	l = size(histEnergy);
	NG = l(1);
	NT = l(2);
	% take the Fourier transform of the electric field on the grid
	NFFT = 2^nextpow2(NG); % Next power of 2 from length of Eg
	Y = fft(histEnergy, NFFT)/NG;
	histSpectrum = 2*abs(Y(1:NFFT/2+1, :));
	histMode = histSpectrum(mode, :);
end