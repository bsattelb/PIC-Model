function growthrate = calc_Two_Stream_Basic(time, histEnergy)
    % Compute critical points and damping rate of energies
    location = find(time == 20);
    log(histEnergy(location));
    growthrate = (log(histEnergy(location)) - log(histEnergy(location-2)))/(time(location) - time(location - 2));
end

