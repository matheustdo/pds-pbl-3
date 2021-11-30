clc;
clear;
close all;

%% Leitura do arquivo de audio
[x, Fs] = audioread('teste1.wav');
x = x(:,1);

%% Determinacao dos parametros do sinal
N = length(x);
Ts = 1/Fs;
dt = N*Ts;
t = 0:Ts:dt-Ts;

%% Geracao do sinal desejado
[d,d_h]= WienerNoiseReduction(x,Fs,1000);


%% Calcula SNR
ruido = x - d;
SNR_TSNR = mag2db(rssq(d(:))/rssq(ruido(:)));

ruido = x - d_h;
SNR_HNRN = mag2db(rssq(d_h(:))/rssq(ruido(:)));
%% Aplicacao dos filtros
ordem = 20;
[y_lms, e_LMS, w_LMS] = LMS1(x, d_h, t, ordem);
[y_rls, e_RLS, w_RLS] = RLS1(x, d_h, t, ordem);

%% Comparação dos erros
figure('Name','Comparação do Erro entre LMS e RLS','NumberTitle','off');
plot(e_LMS);
hold on;
plot(e_RLS);
legend('Erro do LMS','Erro do RLS');
title('Performance de Convergência do LMS e do RLS');

%% Comparação dos sinais
figure('Name','Comparação dos Sinais gerados entre LMS e RLS','NumberTitle','off');
plot(y_rls);
hold on;
plot(y_lms);
legend('Sinal obtido pelo RLS','Sinal obtido pelo LMS');
title('Comparação dos Sinais gerados entre LMS e RLS');

%% Pesos dos filtros LMS e RLS
% figure('Name','Comparação do Pesos entre LMS e RLS','NumberTitle','off');
% stem(w_LMS');
% hold on;
% stem(w_RLS');
% legend('Pesos do LMS','Pesos do RLS');
% title('Comparação do Pesos entre LMS e RLS');

%% Desvio padrão dos erros gerados pelo LMS e RLS
dp_LMS = std(e_LMS);
dp_RLS = std(e_RLS);

%% Geracao dos arquivos de saida
audiowrite('audios/wiener-lms.wav',y_lms,Fs);
audiowrite('audios/wiener-rls.wav',y_rls,Fs);