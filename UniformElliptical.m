close all
clear
clc
A=4; 
B=1;
Pe=[0;0];
Pc=[1; 1];
ang=pi/4;
N=100;
res=zeros(N,1);
dhat=linspace(0,10,N);
M=[cos(ang), -sin(ang);  sin(ang),cos(ang)];
cs=M*(Pc-Pe)+[0;B];
for k=1:N
% syms D1(x)
D1=@(x) cs(2)+sqrt(dhat(k).^2/4-(x-cs(1)).^2);
% syms D2(x)
D2=@(x) cs(2)-sqrt(dhat(k).^2/4-(x-cs(1)).^2);
% syms E1(x)
E1=@(x) B+sqrt(1-(x/A).^2);
% syms E2(x)
E2=@(x) B-sqrt(1-(x/A).^2);


U=@(x) max(E2(x),min(D1(x),E1(x)));
L=@(x) min(E1(x),max(D2(x),E2(x)));
I=@(x) U(x)-L(x);



res(k)=integral(I, max(-A,cs(1)-dhat(k)/2), min(A,cs(1)+dhat(k)/2) );
end
plot(res/A/B/pi)

xlabel('$\hat{d}$', 'Interpreter','Latex')
ylabel('$F_K(\hat{d}\,)$', 'Interpreter','Latex')
% title('CDF for uniform distribution with elliptic support')