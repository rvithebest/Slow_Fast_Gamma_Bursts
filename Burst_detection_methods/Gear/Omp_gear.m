function [b,a,YFIT,res] = Omp_gear(signal,dictSize,sRate,iterMax)
    % make dictionary %
    t_start=tic;
    L=length(signal);
    [Dcos,Dsin,t,f,s] = makeDictionary_mage(dictSize,sRate,signal);
    Y = signal;
    res = [];
    %% OMP starts from here %
    %%%%% Initializations %%
    p = size(Dcos,2);
    J=1:p;
    b = zeros(iterMax,7);
    a = [sRate,L,0,0,0,0,iterMax,0];
    YFIT = zeros(L,1);
    tmpIOPT = zeros(1,iterMax);
    R = signal;
    Xc=Dcos;
    Xs = Dsin;
    %Xc = normc(Xc);
    %Xs=normc(Xs);
    P=[];
    sigE = sum(abs(signal).^2);

    %%%% OMP Loop %%%%

    for k=1:iterMax
        scalProdcos = (R'*Xc);
        scalProdsin = (R'*Xs);
        sqpsum = scalProdsin.^2+scalProdcos.^2;
        [~,i] = max(sqpsum);
        kopt = J(i);
        %J=setdiff(J,kopt);
        %b(k,1) = s(kopt);
        %b(k,2) = f(kopt);
        %b(k,3) = t(kopt);
        my = scalProdsin(i);
        mx = scalProdcos(i);
        phase = atan2(my,mx);
        [gparams] = reassign_gear([t(kopt),f(kopt),2*log(s(kopt))],R,sRate);
        tnw=gparams(1);
        fnw=gparams(2);
        snw=exp(gparams(3)/2);
        atomnew1 = mageAtom_real([tnw,fnw,2*log(snw)],sRate,L,0);
        atomnew2 = mageAtom_imag([tnw,fnw,2*log(snw)],sRate,L,0);
        atomnew1 = normc(atomnew1');
        atomnew2 = normc(atomnew2');
        my1=atomnew2'*R;
        mx1=atomnew1'*R;
        phase1 = atan2(my1,mx1);
        atomnew = mageAtom_real([tnw,fnw,2*log(snw)],sRate,L,-phase1);
        %atomnew = mageAtom_real([t(kopt),f(kopt),2*log(s(kopt))],sRate,L,-phase);
        atomnew = normc(atomnew');
        P=[P atomnew];
        TMP= P\Y;
        R=Y-P*TMP;
        YFIT = P*TMP;
        YFITe = sum(abs(YFIT).^2);
        b(k,1) = snw;
        b(k,2) = fnw;
        b(k,3) = tnw;
        b(k,5) = 1;
        b(k,6) =  -phase1;
        %disp("s0 is "+num2str(s(kopt))+"snw is "+num2str(snw))
        %disp("f0 is "+num2str(f(kopt))+"fnw is "+num2str(fnw))
        %disp("t0 is "+num2str(t(kopt))+"tnw is "+num2str(tnw))
        disp("Reconstruction Energy iteration- "+num2str(k)+" "+num2str((YFITe/sigE)*100))
            if YFITe >= sigE
                break;
            end
            res = [res (YFITe/sigE)*100];
%             clf;
%             plot(signal,'g'); hold on;
%             %plot(R,'b')
%             plot(scalProdcos(i)*normc(mageAtom_real([t(kopt),f(kopt),2*log(s(kopt))],sRate,L,0)'),'k');
%             plot(TMP(k)*atomnew,'r');
%             disp(k);
%             pause;
%         Xc = Dcos(:,J);
%         Xs = Dsin(:,J);
    end
    
    COEFF = P\Y;
    YFIT=P*COEFF;
    for m=1:length(COEFF)
          b(m,4)=COEFF(m);
     end

toc(t_start);
    


end