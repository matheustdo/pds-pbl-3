function LMS1(x, d, t, mu, M)
    N=numel(d);
    wi=zeros(1,M);
    e=[];
    for i=M:N
        e(i)=d(i)-wi*x(i:-1:i-M+1)';
        wi = wi+2*mu*e(i)*x(i:-1:i-M+1);
    end
    y=zeros(N,1);
    for i=M:N
        j=x(i:-1:i-M+1);
        y(i)=((wi)*(j)');
    end
    subplot(221),plot(t,d),title('Sinal desejado'),
    subplot(222),plot(t,x),title('Sinal de entrada ruidoso'),
    subplot(223),plot(t,e),title('Erro'),
    subplot(224),plot(t,y),title('Sinal obtido');

end