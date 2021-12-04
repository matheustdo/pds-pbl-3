function [ha,e,w] = RLS1(x, d, t, M)
    delta = 0.001;
    lambda = 1;
    N = numel(d);
    
    ha = zeros(N,1);
    e = zeros(N,1);
    w = zeros(M,1);
    u = zeros(M,1);
    P = ((1/delta)*eye(M,M));
    for i=1:N
        u = [d(i) 
        u(1:(M-1))];
        x_n = x(i);
        % Step 1: ganho
        g = (((1/lambda)*P*u)/(1+((1/lambda)*u'*P*u)));
        % Step 2: filtragem
        y = (w'*u);
        ha(i) = (ha(i)+y);
        % Step 3: estimando o erro
        E = (x_n-y);
        e(i) = e(i)+E;
        % Step 4: adaptação
        w = w+g*conj(E);
        %Step 5: atualização da correlação
        P = (((1/(lambda))*P)-((1/lambda)*g*u'*P));
    end
    
    textoT2 = strcat('Algoritmo RLS');
    figure('Name',textoT2,'NumberTitle','off');
    subplot(221),plot(t,d),title('Sinal desejado');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    axis([0 t(length(t)) -inf inf]);
    
    subplot(222),plot(t,x),title('Sinal de entrada ruidoso');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    axis([0 t(length(t)) -inf inf]);
    
    subplot(223),plot(t,e),title('Erro');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    axis([0 t(length(t)) -inf inf]);
    
    subplot(224),plot(t,ha),title('Sinal obtido');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    axis([0 t(length(t)) -inf inf]);
end