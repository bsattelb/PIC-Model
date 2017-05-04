function growthrate = calc_Two_Stream_Basic(time, histEnergy, dx)

	histQOI = Energy_Norms.L2Norm(histEnergy, dx);
    % Compute critical points and damping rate of energies
    lnE      = log(histQOI);
    lnE_big  = [lnE(1), lnE, lnE(end)];
    peaks    = find(diff(lnE_big(1:end-1)) > 0 & diff(lnE_big(2:end)) < 0);
    [~, region] = max(diff(histQOI(peaks)));
    startRegion = peaks(region);
    endRegion = peaks(region+1);
    
    cutoff = round((endRegion - startRegion)/10);
    startRegion = startRegion + cutoff;
    endRegion = endRegion - cutoff;

    p = polyfit(time(startRegion:endRegion), lnE(startRegion:endRegion), 1);
    if(p(1) > 0.5)
        semilogy(time(2:end), histQOI, time(2:end), exp(time(2:end)*p(1) + p(2)));
        p(1)
        exit()
    end
    growthrate = p(1);
end

