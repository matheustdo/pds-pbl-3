clc;
clear;
close all;
N = 10000;
t = 0:N-1;
d = sin([1:N]*0.05*pi); % sinal desejado
ruido  = randn(1,N)*0.1; % ruído 
v1= filter(1,[1 -0.8],ruido); % ruído filtrado
x = d + v1;
mu = 0.002;
ordem = 20;
valorSNR = mag2db(rssq(d(:))/rssq(ruido(:)));
n0 = 10; % delay de 10 amostras
len = N - n0; % tamanho do vetor reduzido
x_del = zeros(N,1);

for i = 1:len
    x_del(i) = x(i+n0);
end

% figure(1);
% plot(x(1:1000),':k')
% hold on
% plot(x_del(1:1000),'k')
% legend('x(n)', 'x(n-n0)')
% title('Reference signal x(n-n0)');
% xlabel('samples n')
% ylabel('amplitude')
% grid on

% LMS1(x, x_del, t, ordem);
% RLS1(x, x_del, t, ordem);
LMS1(x, d, t, ordem);
RLS1(x, d, t, ordem);
ha = adaptfilt.lms(ordem,mu);
[y,e_M] = filter(ha,x,x_del);
valorSNR_M = snr(d,e_M);
textoT = strcat('Algoritmo LMS do MATLAB | SNR = ',num2str(valorSNR_M));
figure('Name',textoT,'NumberTitle','off');
subplot(221),plot(t,x_del),title('Sinal desejado') , axis([0 1000 -3 3]);
subplot(222),plot(t,x),title('Sinal de entrada ruidoso'), axis([0 1000 -3 3]);
subplot(223),plot(t,e_M),title('Erro'), axis([0 1000 -3 3]);
subplot(224),plot(t,y),title('Sinal obtido'); axis([0 1000 -3 3]);

ffR = 1;
rls = adaptfilt.rls(ordem,ffR);
[y_RLS, e_RLS] = filter(rls,x,x_del);
textoT2 = strcat('Algoritmo RLS do MATLAB | SNR = ',num2str(valorSNR));
figure('Name',textoT2,'NumberTitle','off');
subplot(221),plot(t,x_del),title('Sinal desejado'),axis([0 1000 -3 3]);
subplot(222),plot(t,x),title('Sinal de entrada ruidoso'),axis([0 1000 -3 3]);
subplot(223),plot(t,e_RLS),title('Erro'),axis([0 1000 -3 3]);
subplot(224),plot(t,y_RLS),title('Sinal obtido');axis([0 1000 -3 3]);