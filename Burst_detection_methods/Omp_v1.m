function [batoms,header,YFIT] = doOMPv2(analogData,timeVals,trial,N1,itermax)
L = length(analogData); 
sRate = 1/(timeVals(2)-timeVals(1));
%% Dictionary Generation %%            
N = N1; 
Dcos = zeros(L,N1);
Dsin = zeros(L,N1);
rng shuffle
s = randi([1 L],[1 N]);
s1 = s/sRate;
p1 = randperm(N,N1);
s = s(p1);
s1 = s1(p1);
clear p1;
rng shuffle
f = randi([1 (L/2)],[1 N]);
f1 = f*sRate/L;
p2 = randperm(N,N1);
f1 = f1(p2);
f=f(p2);
clear p2;
rng shuffle
t = randi([0 L],[1 N]);
t1 = t/sRate;
p3 = randperm(N,N1);
t=t(p3);
t1=t1(p3);
clear p3;
clear s;
clear t;
clear f;
for i=1:N1
    Dcos(:,i) = gaborCos(s1(i),f1(i),t1(i),timeVals,0);
    Dsin(:,i) = gaborSin(s1(i),f1(i),t1(i),timeVals,0);
end

% Normalization of columns of Dcos an Dsin.
S1 = sum(Dcos.*Dcos).^0.5;
Dcos = Dcos./repmat(S1,L,1);% the columns norm are set to 1 
S2 = sum(Dsin.*Dsin).^0.5;
Dsin = Dsin./repmat(S2,L,1);% the columns norm are set to 1
Y = analogData;

            


%% OMP starts frome here %%
p = size(Dcos,2); % number of atoms in the dictionary
J     = 1:p;                                  % index of remaining vectors
batoms = zeros(itermax,7);% in the dictionary
header = [sRate,L,0,0,0,0,itermax,0];
YFIT  = zeros(L,1) ;
%%tmpqual  = zeros(1,itermax);
tmpIOPT  = zeros(1,itermax);                  % Index of the selected vectors
R = analogData;
%P = zeros(1024,itermax);
Xc  = Dcos;
Xs = Dsin;
%ss =s1;
%ff =f1;
%tt =t1;
P=[];

        for k = 1:itermax
            scalProdcos = (R' * Xc);
            scalProdsin  = (R' * Xs);
            sqpsum = scalProdsin.^2+scalProdcos.^2;
            % OMP and ...
            [~, i] = max(sqpsum);
            kopt = J(i);                           % index of the kept atom.
            J    = setdiff(J,kopt);
            tmpIOPT(k)	= kopt;                    % update support
            batoms(k,1) = s1(kopt);
            batoms(k,2) = f1(kopt);
            batoms(k,3) = t1(kopt);
            % phase fincding step %
            my = scalProdsin(i);
            mx = scalProdcos(i);
            phase = atan2(my,mx);
            batoms(k,5)=phase;
            atomnew = gaborCos(s1(kopt),f1(kopt),t1(kopt),timeVals,-phase);
            atomnew = atomnew/(sum(atomnew.*atomnew).^0.5);
            P= [P atomnew'];
            TMP  = P\Y;                            % least square estimate
            R = Y - P*TMP;
            YFIT = P*TMP;  
            
            Xc = Dcos(:,J);
            Xs = Dsin(:,J);
           % ss = s1(J);
           % ff = f1(J);
           % tt = t1(J);
            
        end

        COEFF = P\Y;
        YFIT  = P * COEFF;
        batoms(:,4) = TMP;
        %disp(norm(R)/norm(Y));
end
function array = gaborCos(s,f,u,timeVals,phase)
    signal_sampling=(timeVals(2)-timeVals(1));
    signal_size = length(timeVals);
    t = (0:signal_size-1)*signal_sampling;
    array = exp(-pi*((t-u)/s).^2).*cos(2*pi*f.*(t)+phase);
end       
function array = gaborSin(s,f,u,timeVals,phase)
    signal_sampling=(timeVals(2)-timeVals(1));
    signal_size = length(timeVals);
    t = (0:signal_size-1)*signal_sampling;
    array = exp(-pi*((t-u)/s).^2).*sin(2*pi*f.*(t)+phase);
end      
