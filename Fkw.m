clear
clc
close all

%Numerical integration problems for lots of cases - integral2 yeilds NaN.
% Pbar=rand([2,1])*8-5;
% a=Pbar(1);
% b=Pbar(2);
% R=sqrt(a^2+b^2);
% sx=rand(1);
% sy=rand(1);
% rho=rand(1);

a=1;
b=1;
Pbar=[a;b];
R=sqrt(a^2+b^2);
sx=1;
sy=1/2;
rho=0;                          %==0 works with different sx sy; -10 works; -7 gives NaN on first only
Sigma=[sx^2 rho*sx*sy;rho*sx*sy, sy^2];     
% Sigma=3*eye(2);
syms fxy(x,y)
fxy(x,y)=mvnpdf([x;y], Pbar,Sigma);
fsurf(fxy(x,y))
syms gamma1(z,w)
gamma1(z,w)=(1/2-(z-w+1)/(2*R^2));
syms gamma2(z,w)
gamma2(z,w)=1/2*sqrt(2*(z+w-1)/R^2-(z-w+1).^2/R^4-1);

syms x1(z,w)
x1(z,w)=gamma1(z,w)*a-gamma2(z,w)*b;
syms y1(z,w)
y1(z,w)=gamma1(z,w)*b+gamma2(z,w)*a;
syms x2(z,w)
x2(z,w)=gamma1(z,w)*a+gamma2(z,w)*b;
syms y2(z,w)
y2(z,w)=gamma1(z,w)*b-gamma2(z,w)*a;


syms fzw(z,w)
fzw(z,w)=1/(4*R^2*gamma2(z,w))*(eval(fxy(x1(z,w),y1(z,w)))+eval(fxy(x2(z,w),y2(z,w))));
f=@(w,z) eval(fzw(z,w));            %swap order because integral2 integrates y first


syms fzw_ext(w,z)
fzw_ext(z,w)=piecewise(abs(R-sqrt(z))<sqrt(w-1),fzw(x,y),0);
f2=@(w,z) eval(fzw_ext(z,w)); 
%% Optional section - verification that complete integral is 1.
lb=@(x) (R-sqrt(x-1)).^2;
ub=@(x) (R+sqrt(x-1)).^2;
integral2(f,1, inf,lb,ub, 'Method','iterated')

%% f2 trial
% assume(w-1, 'positive')
% assume(z, 'positive')
% n=5;
% Fkw=zeros(1,n);
% ds=linspace(0,1,n);
% for k=1:n
%     d=ds(k)
%     ub=@(x) x*d.^2.*(1+R.^2);
%     Fkw(k)=integral2(f2,1, inf,0,ub, 'Method','iterated')
% %     subplot(n/2,2, k)
% %     plot(t,(R-sqrt(t-1)).^2, 'b')
% %     hold on
% %     plot(t, t*d.^2.*(1+R.^2), 'k')
% %     plot(t,(R+sqrt(t-1)).^2, 'r')
% %     subtitle(sprintf('d=%f',d))
%     
% end
% % legend('(R-\sqrt(t-1)).^2','t*d.^2.*(1+R.^2)','(R+\sqrt(t-1)).^2')
% Fkw(isnan(Fkw))=0;


%% OBS takes quite a while to run even for n=25

n=25;
Fkw=zeros(1,n);
ds=linspace(0,1,n);
t=linspace(1,100,100);


for k=1:n
    d=ds(k)
    lb=@(x) min(x*d.^2.*(1+R.^2),(R-sqrt(x-1)).^2);     %if xd^2(1+R^2)<(R-sqrt(x-1))^2<(R+sqrt(x-1))^2 then integral 
                                                        %becomes from xd^2(1+R^2) to xd^2(1+R^2) (should be 0)
    ub=@(x) min(x*d.^2.*(1+R.^2),(R+sqrt(x-1)).^2);
    Fkw(k)=integral2(f,1, inf,lb,ub, 'Method','iterated')
%     subplot(n/2,2, k)
%     plot(t,(R-sqrt(t-1)).^2, 'b')
%     hold on
%     plot(t, t*d.^2.*(1+R.^2), 'k')
%     plot(t,(R+sqrt(t-1)).^2, 'r')
%     subtitle(sprintf('d=%f',d))
%     if isnan(Fkw(k))
%         'NaN'
%         lim=@(x) (R+sqrt(x-1)).^2;
% %         lim2=@(x) (R-sqrt(x-1)).^2;
%         Fkw(k)=1-integral2(f,1, inf,ub,lim, 'Method','iterated')
%     end
    
end
% legend('(R-\sqrt(t-1)).^2','t*d.^2.*(1+R.^2)','(R+\sqrt(t-1)).^2')
Fkw(isnan(Fkw))=0;

%% Plot result
figure(1),
plot(linspace(0,1,n),Fkw);
title('CDF of \kappa_\omega for a gaussian distibution P_G')
xlabel('d')
ylabel('F_{\kappa_\omega} (d)')
ylim([0,1])


