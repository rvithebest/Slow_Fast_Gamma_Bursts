function [dictCos,dictSin,tVals,fVals,sVals] = makeDictionary_mage(dictSize,sRate,signal)
    L = length(signal);
    dictCos = zeros(L,dictSize);
    dictSin = zeros(L,dictSize);
    rng shuffle
    fVals = randi([1 (L/2)],[1 dictSize]);
    fVals = fVals*sRate/L;
    rng shuffle
    tVals = randi([0 L],[1 dictSize]);
    tVals = tVals/sRate;
    rng shuffle
    sVals = randi([0 L-1],[1 dictSize]);
    
    sVals = sVals/250;
    %  dictionary generation %%
    for i=1:dictSize
        f_tmp = fVals(i)+rand*sRate/L;
        fVals(i)=f_tmp;
        t_tmp = tVals(i)+rand/sRate;
        tVals(i)=t_tmp;
        s_tmp = sVals(i)+rand/sRate;
        sVals(i)=s_tmp;
        dictCos(:,i) = mageAtom_real([t_tmp,f_tmp,2*log(s_tmp)],sRate,L,0);
        dictSin(:,i) = mageAtom_imag([t_tmp,f_tmp,2*log(s_tmp)],sRate,L,0);


        
    
    
    end
    %normalize the dictionary %%
    S1 = sum(dictCos.*dictCos).^0.5;
    dictCos = dictCos./repmat(S1,L,1);
    S2 = sum(dictSin.*dictSin).^0.5;
    dictSin = dictSin./repmat(S2,L,1);
    dictCos = normc(dictCos);
    dictSin = normc(dictSin);

end