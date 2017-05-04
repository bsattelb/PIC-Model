function histEL2 = L2Norm(histEnergy, dx)
	% Potential energy
	Efield = 0.5*sum(histEnergy.^2)*dx;
	%L2 norm
	histEL2 = sqrt(2*Efield);
end