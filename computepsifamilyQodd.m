function [win,psi,lpsi,time,cycles]=computepsifamilyQodd(F,sp,cycles1,cycles2,winsize)

% Inputs:
%
% F is a vector of frequencies. A step of 1 Hz is suggested. (e.g.,
% F=1:100).
%
% sp is the sampling step or 1./sampling frequency. (e.g., 1/256).
%
% cycles1 is the number of cycles for the lowest frequency (e.g.,
% cycles1=1);
%
% cycles2 is the number of cycles for the highest frequency (e.g.,
% cycles2=20);
%
% winsize is the number of samples of your wavelet. A simple recomendation
% is to set it as:
% winsize=max(pow2(nextpow2(frames)-4),4)+1; 
% Note: the 'plus one' is to have odd-length wavelets, which allows for perfectly centered convolutions.
% 
%
% Outputs:
%
% win is a matrix with wavelets in each row, so that the wavelet transform
% is computed by W=win*signal; (signal being a column vector with the same sampling period).


lf=length(F);
if cycles1==cycles2
    cycles=cycles1*ones(1,lf);
else
    cycles=cycles1:(cycles2-cycles1)/(length(F)-1):cycles2;
end
win=zeros(lf,winsize);
for k=1:lf
	fk=F(k);
	sigf=fk/cycles(k);
	sigt=1./(2*pi*sigf);
	tneg=[-sp:-sp:-sp*winsize/2];
	tpos=[0:sp:sp*winsize/2];
	t=[fliplr(tneg) tpos];
	A=1./sqrt(sigt*sqrt(pi));
	p=A.*(exp(-(t.^2)./(2*(sigt^2))).*exp(2*i*pi*fk*t));
    psi{k}=p;
    lpsi(k)=length(psi{k});
	time{k}=t;
	win(k,:)=p;
end