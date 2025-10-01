
function array = gaborCos(s,f,u,timeVals,phase)
    signal_sampling=(timeVals(2)-timeVals(1));
    signal_size = length(timeVals);
    t = (0:signal_size-1)*signal_sampling;
    array = exp(-pi*((t-u)/s).^2).*cos(2*pi*f.*(t)+phase);
end       