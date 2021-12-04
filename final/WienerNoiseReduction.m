function [esTSNR,esHRNR]=WienerNoiseReduction(ns,fs,IS)

% [esTSNR,esHRNR]=WIENERNOISEREDUCTION(ns,fs,IS) 
%
% Parametros de entrada :  
%   ns          Sinal ruidoso
%   fs          Frequencia de amostragem (in Hz)
%   IS          número de amostras do período inicial de atividade sem fala ou IS (ou seja, período de silêncio inicial)  
%               
%
% Parametros de saída :
%   esTSNR      fala aprimorada com o método de redução de ruído em duas
%   etapas (TSNR)
%
%   esHNRN      discurso aprimorado com o método de redução de ruído de regeneração harmônica 

%

%%  Entrada Ruidosa 

l = length(ns);
s=ns;

wl = fix(0.025*fs);    % o comprimento da janela é de 25 ms 
NFFT=2*wl;             % O tamanho da FFT é o dobro do comprimento da janela 
hanwin = hanning(wl);


if (nargin<3 | isstruct(IS))
    IS=10*wl;             %Silêncio inicial ou apenas parte do ruído nas amostras (= dez quadros) 
end
%%  Computa as estatísticas do ruído


nsum = zeros(NFFT,1);

count = 0; 
    for m = 0:IS-wl
     nwin = ns(m+1:m+wl).*hanwin;	
      nsum = nsum + abs(fft(nwin,NFFT)).^2;
     count = count + 1;
    end

d= (nsum)/count;


%%  Parte Principal
SP = 0.4;      % A porcentagem de deslocamento é de 60%. O método de sobreposição-adição funciona bem com este valor 
normFactor=1/SP;
overlap = fix((1-SP)*wl); % sobreposição entre quadros sucessivos 
offset = wl - overlap;
max_m = fix((l-NFFT)/offset);

zvector = zeros(NFFT,1);
oldmag = zeros(NFFT,1);
news = zeros(l,1);

phasea=zeros(NFFT,max_m);
xmaga=zeros(NFFT,max_m);
tsnra=zeros(NFFT,max_m);
newmags=zeros(NFFT,max_m);

alpha = 0.99;

%Iteração para remover ruído 

%%  TSNR
for m = 0:max_m
   begin = m*offset+1;    
   iend = m*offset+wl;
   speech = ns(begin:iend);       %Extração da parte da fala
   winy = hanwin.*speech;   %Executa a janela de hanning
   ffty = fft(winy,NFFT);          %Realiza a FFT
   phasey = angle(ffty);         %Pega a fase do sinal
   phasea(:,m+1)=phasey;       %para uso do HRNR 
   magy = abs(ffty);             %Pega a magnitude do sinal
   xmaga(:,m+1)= magy;           %para uso do HRNR 
   postsnr = ((magy.^2) ./ d)-1 ;      %calcula a posteriori SNR 
   postsnr=max(postsnr,0.1);  % limitação para evitar distorção 
   
   %calcula a priori SNR usando abordagem direcionada à decisão 
   eta = alpha * ( (oldmag.^2)./d ) + (1-alpha) * postsnr;
   newmag = (eta./(eta+1)).*  magy;
   
   %Calcula TSNR
   tsnr = (newmag.^2) ./ d;
   Gtsnr = tsnr ./ (tsnr+1);         %ganho do TSNR 
   tsnra(:,m+1)=Gtsnr;    
   Gtsnr=max(Gtsnr,0.15);  
   Gtsnr = gaincontrol(Gtsnr,NFFT/2);
   
      %para uso do HRNR
   newmag = Gtsnr .* magy;
   newmags(:,m+1) = newmag;     %para uso do HRNR
   ffty = newmag.*exp(i*phasey);
   oldmag = abs(newmag);
   news(begin:begin+NFFT-1) = news(begin:begin+NFFT-1) + real(ifft(ffty,NFFT))/normFactor;
end

esTSNR=news;


%%  HRNR
%Não-linearidade
newharm= max(esTSNR,0);
news = zeros(l,1);
%
for m = 0:max_m
   begin = m*offset+1;    
   iend = m*offset+wl;

   nharm = hanwin.*newharm(begin:iend);
   ffth = abs(fft(nharm,NFFT));          %Executa uma FFT

   snrham= ( (tsnra(:,m+1)).*(abs(newmags(:,m+1)).^2) + (1-(tsnra(:,m+1))) .* (ffth.^2) ) ./d;
   
   newgain= (snrham./(snrham+1));
   newgain=max(newgain,0.15);  
   
   newgain=gaincontrol(newgain,NFFT/2);
   
   newmag = newgain .*  xmaga(:,m+1);
 
   ffty = newmag.*exp(i*phasea(:,m+1));
   
   news(begin:begin+NFFT-1) = news(begin:begin+NFFT-1) + real(ifft(ffty,NFFT))/normFactor;
end;

%Saida
esHRNR=news;
esTSNR = esTSNR * max(abs(ns))/max(abs(esTSNR));
esHRNR = esHRNR * max(abs(ns))/max(abs(esHRNR));

figure;
[B,f,T] = specgram(ns,NFFT,fs,hanning(wl),wl-10);
imagesc(T,f,20*log10(abs(B)));axis xy;colorbar
title(['Espectro do Sinal - Sinal Ruidoso']);
xlabel('Tempo (segundos)');ylabel('Frequencia (Hz)');

figure;
[B,f,T] = specgram(esTSNR,NFFT,fs,hanning(wl),wl-10);
imagesc(T,f,20*log10(abs(B)));axis xy;colorbar
title(['Espectro do Sinal - Saída de voz TSNR ']);
xlabel('Tempo (segundos)');ylabel('Frequencia (Hz)');

figure;
[B,f,T] = specgram(esHRNR,NFFT,fs,hanning(wl),wl-10);
imagesc(T,f,20*log10(abs(B)));axis xy;colorbar
title(['Espectro do Sinal - Saída de voz HRNR']);
xlabel('Tempo (segundos)');ylabel('Frequencia (Hz)');




function        NewGain=gaincontrol(Gain,ConstraintInLength)
%
%Title  : Additional Constraint on the impulse response  
%         to ensure linear convolution property
%
%
%Description : 
%
% 1- The time-duration of noisy speech frame is equal to L1 samples.
%
% 2- This frame is then converted in the frequency domain 
%       by applying a short-time Fourier transform of size NFFT leading
%       to X(wk) k=0,...,NFFT-1 when NFFT is the FFT size.
%
% 3- The estimated noise reduction filter is G(wk) k=0,1,...,NFFT-1 
%       leading to an equivalent impulse response g(n)=IFFT[G(wk)] 
%       of length L2=NFFT
%
% 4- When applying the noise reduction filter G(wk) to the noisy 
%       speech spectrum X(wk), the multiplication S(wk)=G(wk)X(wk) is
%       equivalent to a convolution in the time domain. So the
%       time-duration of the enhanced speech s(n) should be equal to 
%       Ltot=L1+L2-1.
%
% 5- If the length Ltot is greater than the time-duration of the IFFT[S(wk)] 
%       the a time-aliasing effect will appear.
%
% 6- To overcome this phenomenon, the time-duration L2 of the equivalent
%       impulse response g(n) should be chosen such that Ltot = L1 + L2 -1 <= NFFT 
%       => L2 <= NFFT+1-Ll
%
%       here we have NFFT=2*Ll so we should have L2 <= Ll+1. I have made
%       the following choice : the time-duration of g(n) is limited to
%       L2=NFFT/2=L1 (see lines 88 and 192)
%
%Author : SCALART Pascal
%
%October  2008
%


meanGain=mean(Gain.^2);
NFFT=length(Gain);

L2=ConstraintInLength;

win=hamming(L2);

% Frequency -> Time
% computation of the non-constrained impulse response
ImpulseR=real(ifft(Gain));

% application of the constraint in the time domain
ImpulseR2=[ImpulseR(1:L2/2).*win(1+L2/2:L2);zeros(NFFT-L2,1);ImpulseR(NFFT-L2/2+1:NFFT).*win(1:L2/2)];

% Time -> Frequency
NewGain=abs(fft(ImpulseR2,NFFT));

meanNewGain=mean(NewGain.^2);

NewGain=NewGain*sqrt(meanGain/meanNewGain); % normalisation to keep the same energy (if white r.v.)


   
