function ...
    gparams_out=mage_reassign_gparams_td(rs,srate,gparams_in)

%  gparams_out=mage_reassign_gparams_td(rs,srate,gparams_in)

% Author: Ryan T. Canolty (2019) rcanolty@gmail.com http://accl.psy.vanderbilt.edu/
% For additional information, visit https://www.biorxiv.org to read:
% TITLE: Multiscale Adaptive Gabor Expansion (MAGE):
% Improved Detection of Transient Oscillatory Burst Amplitude and Phase
% AUTHORS: Ryan T. Canolty and Thilo Womelsdorf  
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


rs=rs(:);
Nt=length(rs);

tp=gparams_in(:,1);
vp=gparams_in(:,2);
sp=gparams_in(:,3);
% tstd=sqrt(exp(sp)/(4*pi));

t=(0:Nt-1)/srate;
t=t(:);

gp=exp(-exp(-sp).*pi*(t-tp).^2).*... % amplitude envelope
        exp(2*pi*1i*vp.*(t-tp));    % linear frequency chirp
gp=gp/norm(gp); % need to normalize due to discrete sampling

pre1=2*exp(-sp)*pi*(t-tp) - 2*pi*1i*vp;
pre2=2*pi*1i*(t-tp);
pre3=(-1/4)+exp(-sp)*pi*(t-tp).^2;

dgpdtp=pre1.*gp;
dgpdvp=pre2.*gp;
dgpdsp=pre3.*gp;

clear fs
fs(1)=sum(rs.*conj(gp));
fs(2)=sum(rs.*conj(dgpdtp));
fs(3)=sum(rs.*conj(dgpdvp));
fs(4)=sum(rs.*conj(dgpdsp));

rztp=real(fs(2)./fs(1))';
rzvp=real(fs(3)./fs(1))';
rzsp=real(fs(4)./fs(1))';

st_est=sp+real(log(-1+(2*pi*exp(sp))./...
    (-rzvp.^2+exp(sp).*(exp(sp).*rztp.^2+pi-4*pi*rzsp))));

vt_est=vp+((exp(-sp)+exp(-st_est))./(2*pi)).*rzvp;

tt_est=tp+((exp(sp)+exp(st_est))./(2*pi)).*rztp;

gt_est=exp(-exp(-st_est)*pi*(t-tt_est).^2).*... % amplitude envelope
        exp(2*pi*1i*vt_est*(t-tt_est));    % linear frequency chirp
gt_est=gt_est/norm(gt_est); % need to normalize due to discrete sampling
val=sum(rs.*conj(gt_est));
% val=1;
gparams_out=[tt_est vt_est st_est 2*abs(val) angle(val)];



