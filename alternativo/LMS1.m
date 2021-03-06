function [ha,e,w] = LMS1(x, d, t, M)
    mu = 0.005;
    N=numel(d);
    ha = zeros(N,1);
    e = zeros(N,1);
    w = zeros(M,1);
    U = zeros(M,1);

    % Step 1: Filtragem
    for i=1:N
        U = [d(i) 
        U(1:(M-1))];
        x_n = x(i);
        y = (w'*U);
        ha(i) = (ha(i)+y);
        % Step 2: Estimando o erro
        E_LMS = (x_n-y);
        e(i) = e(i)+E_LMS;
        % Step 3: Adapta??o
        w = (w+(mu*E_LMS*U));
    end
    ruido = x - d;
    valorSNR = mag2db(rssq(d(:))/rssq(ruido(:)));
    textoT2 = strcat('Algoritmo LMS | SNR = ',num2str(valorSNR));
    figure('Name',textoT2,'NumberTitle','off');
    subplot(221),plot(t,d),title('Sinal desejado');
    subplot(222),plot(t,x),title('Sinal de entrada ruidoso');
    subplot(223),plot(t,e),title('Erro');
    subplot(224),plot(t,ha),title('Sinal obtido');
    %axis([0 1000 -3 3]);
end