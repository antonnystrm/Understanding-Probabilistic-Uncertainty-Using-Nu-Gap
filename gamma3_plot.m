clear
clc

a=2;
b=2;
R=norm([a,b]);
syms gamma3_ext(z,w)
% gamma3_ext(z,w)=piecewise(abs(R-sqrt(z))<sqrt(w-1), sqrt(8*R^2*(z+w-1)-4*(z-w+1).^2-4*R^4), 0);
gamma3_ext(z,w)=piecewise(abs(R-sqrt(z))<sqrt(w-1), 1/2*sqrt(2/R^2*(z+w-1)-(z-w+1).^2/R^4-1), 0);
%gamma2

%% plot gamma3 and boundaries
h=figure(1);
xlim([0 100])
hold on
syms fill(z,w)
fill(z,w)=0;
fsurf(fill, [0 100 1 100],'MeshDensity',20);
% title('$c_2$ plotted on the region $|R-\sqrt{l}|<\sqrt{t-1}$,\; $r=2*\sqrt{2}$', 'Interpreter', 'Latex')xlabel('l',"Interpreter", "Latex")
xlabel('$l$',"Interpreter", "Latex")
ylabel('$t$',"Interpreter", "Latex")
y=linspace(1,100,100);
plot((R-sqrt(y-1)).^2,y,'b','LineWidth',5)
plot((R+sqrt(y-1)).^2,y,'b','LineWidth',5)
% plot([0 R^2/4, R^2], [R^2+1, R^2/4+1 1], 'ro','LineWidth',5)
xlim([0,100])
fsurf(gamma3_ext, [0 100 1 100],'MeshDensity',20);
plot3(y-1,y,3.5*ones(length(y),1),'r','LineWidth',5);
colorbar

% z=linspace(0,10,100);
% w=linspace(1,11,100);
% [Z,W]=meshgrid(linspace(0,10,5),linspace(1,11,5));

%figure(1)
%surf(Z,W,abs(gamma3_ext(Z,W)))
hold off


%% form and plot fzw
Pbar=[a;b];
Sigma=eye(2);
syms fxy(x,y)
fxy(x,y)=mvnpdf([x,y], Pbar',Sigma);
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
fzw(z,w)=1/gamma3_ext(z,w)*(eval(fxy(x1(z,w),y1(z,w)))+eval(fxy(x2(z,w),y2(z,w))));
eval(fzw(1,2))  
% figure(2)
% fsurf(fzw)%, [0 10 1 11])         BROKEN for some reason see next section

% 


% syms x
% y = piecewise(x < -2,-2,-2 < x < 2,x,x > 2,2);
% fplot(y)
% syms f(x,y)
% f(x,y) = real(atan(x + i*y));
% fsurf(f)
%%
% z=linspace(0,10,100);
% w=linspace(1,11,100);
[Z,W]=meshgrid(linspace(0,5,50),linspace(1,6,50));

figure(3)
surf(Z,W,real(eval(fzw(Z,W))))
title('Joint PDF f_{zw}(z,w)')
xlabel('z')
ylabel('w')
hold on
y=linspace(1,6,100);
plot((R-sqrt(y-1)).^2,y,'b','LineWidth',5)
plot((R+sqrt(y-1)).^2,y,'b','LineWidth',5)
xlim([0 5])
zlim([0 1])
hold off

%% integral     BROKEN
% Fzw=int(fzw);
% figure(4)
% surf(Z,W,real(eval(Fzw(Z,W))))



%%
%plot whats under root to see when negative.
figure(5)
[Z,W]=meshgrid(linspace(-1,10,100),linspace(-1,10,100));
surf(Z,W,8*R^2*(Z+W-1)-4*(Z-W+1).^2-4*R^4)
hold on
surf(Z,W,ones(size(Z)))
y=linspace(1,10,100);
plot((R-sqrt(y-1)).^2,y,'b','LineWidth',5)
plot((R+sqrt(y-1)).^2,y,'b','LineWidth',5)
xlim([-1 10])
hold off
title('inside root of \gamma_3')
xlabel('z')
ylabel('w')
%negative -> complex so shouldnt be an issue.

