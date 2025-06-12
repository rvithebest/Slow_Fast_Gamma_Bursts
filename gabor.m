%%% This function is modified from the original version provided by Prof Piotr durka
%%% Subhash Chandran


function s=gabor(signal_size, signal_sampling, width, frequency, position, amplitude, phase)
%s=gabor(signal_size, signal_sampling, width, frequency, position, amplitude, phase)
% all values in physical units: seconds, Hz etc.

t=0:signal_size*signal_sampling-1;
%parameters in samples (points):
p=position*signal_sampling;
w=width*signal_sampling;
f=frequency/signal_sampling;

if width==signal_size
    s=(amplitude/2)*cos(2*pi*f.*t+phase);
elseif width==0
    s=zeros(size(t));s(round(p)+1)=(amplitude/2);
elseif frequency==0
    s=(amplitude/2)*exp(-pi*((t-p)/w).^2);
else
    s=(amplitude/2)*exp(-pi*((t-p)/w).^2).*cos(2*pi*f.*t+phase);
end

%...
%subplot 211;plot(t./signal_sampling,s);
%subplot 212; spectrum(s,[],[],[],signal_sampling); set(gca,'yscale','linear')
