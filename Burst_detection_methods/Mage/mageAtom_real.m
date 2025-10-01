function atom = mageAtom_real(params,sRate,L,phase)
    t = params(1);
    v=params(2);
    s=params(3);
    timeVals = (0:L-1)*(1/sRate);
    atom = ((2)^(1/4)).*exp(-(s/4)-(pi*exp(-s)*(timeVals-t).^2)).*cos((2*pi*v*(timeVals-t))+phase);
end