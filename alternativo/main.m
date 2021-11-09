clc;
clear;
close all;
N = 2000;
t = 0:N-1;
w0 = 0.01;  
d = sin(2*pi*[1:N]*w0);
ruido = randn(1,N)*0.1;
x = d + ruido;
mu = 0.005;
ordem = 30;
valorSNR = mag2db(rssq(d(:))/rssq(ruido(:)));
% valorSNR = snr(d,ruido);
LMS1(x,d,t,mu,ordem,valorSNR);
ha = adaptfilt.lms(ordem,mu);
[y,e_M] = filter(ha,x,d);
textoT = strcat('Algoritmo LMS do MATLAB | SNR = ',num2str(valorSNR));
figure('Name',textoT,'NumberTitle','off');
subplot(221),plot(t,d),title('Sinal desejado'),
subplot(222),plot(t,x),title('Sinal de entrada ruidoso'),
subplot(223),plot(t,e_M),title('Erro'),
subplot(224),plot(t,y),title('Sinal obtido');




