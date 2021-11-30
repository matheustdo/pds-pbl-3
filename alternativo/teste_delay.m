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
n_delay = 10; % numero de amostras do atraso
d_len = N - n_delay; % tamanho do vetor reduzido
d = zeros(N,1);
for i = 1:d_len
    d(i) = x(i+n_delay);
end

%% Aplicacao dos filtros
ordem = 20;
y_lms = LMS1(x, d, t, ordem);
y_rls = RLS1(x, d, t, ordem);

%% Geracao dos arquivos de saida
audiowrite('audios/delay-lms.wav',y_lms,Fs);
audiowrite('audios/delay-rls.wav',y_rls,Fs);
