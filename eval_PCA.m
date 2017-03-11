clear 
close all

[x,Fs]=wavread('Opptak_1.wav');
n=length(x)
figure
plot((1:length(x))/Fs,x)
%resampler for å simulere lavere samplingsfrekvens
Fs_old=Fs;
Fs=20000;
x=resample(x,Fs,Fs_old);

%log-spectrogram, benytter blacmannvindu med vinduslengde 128
n_hopp=16%antall samples som vindu forskyves med for hver SFT
for kk=1:((floor(length(x)-128+n_hopp)/n_hopp))
    temp=log10(abs(fft(blackman(128).*x((1:128)+(kk-1)*n_hopp))));
    %utelater speilspektrum
S(:,kk)=temp(1:end/2+1,:);
end
S2=S(:,end/2+1:end);
S1=S(:,1:end/2);
[~,n]=size(S);
[~,n2]=size(S1);
% estimerer kanal transferfunksjon for første del av tidssekvensen
S_channel=mean(S1,2);
%kompenserer for estimert kanal transferfunksjon
S1=S1-S_channel(:,ones(length(S1),1));
S2=S2-S_channel(:,ones(length(S2),1));
S=S-S_channel(:,ones(length(S),1));

 
 figure
pcolor((1:n)*32/Fs,(0:64)*Fs/128,S);
colorbar
caxis([0 2])
shading flat
title('log-spectrogram')
F=(0:128/2)*Fs/128;
figure
title('channel log-spectrum response estimat')
plot(F,S_channel);
%PCA basert på første del av tidssekvensen
CC=S1*S1';
 [U,S_,V] = svd((CC));

  figure
pcolor((1:n)*32/Fs,1:length(F),V'*S)
colorbar
caxis([-2 2])
title('PCA')
shading flat

 eigen_val=diag(S_);
 figure
 plot((eigen_val))
 V=V(:,1:10);

 figure
pcolor((1:n)*32/Fs,1:10,V'*S)
colorbar
caxis([-2 2])
title('PCA')
shading flat

  figure
pcolor(((n-n2+1):n)*32/Fs,1:10,V'*S2)
colorbar
caxis([-2 2])
title('PCA the second part of time sequence')
shading flat
