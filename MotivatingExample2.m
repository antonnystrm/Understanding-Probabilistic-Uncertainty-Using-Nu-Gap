clear
clc
s=tf('s');
G_nom=2/(1+s/10)^3;         %Nominal model

%generate some io data
h1 = 0.1;
Gd = c2d(G_nom,h1)
N = 1000;
frequencies = [0 0.05];
u = idinput(N,'rbs',frequencies);
[y,t] = lsim(Gd,u);
noiselevel = 0.1*std(y);
e = noiselevel*randn(N,1);
yn = y + e;

noisydata = iddata(yn,u,h1);
% noisefreedata = iddata(y,u,h1);


%identify with correct model structure in continuous time
nz=0;np=3;
G_est=tfest(noisydata,np, nz)
%resid(noisefreedata,G_est)
%%
figure(1),
w=linspace(0,99,1000);
h1=nyquistplot(G_est,w,'k');
setoptions(h1,'ConfidenceRegionDisplaySpacing',1,'ShowFullContour','off','ConfidenceRegionNumberSD',10);
hold on
plot(-1,0,'r+','MarkerSize',15)
% zoomcp(h)
showConfidence(h1)
hold off
%optional rerun with new color.

setoptions(h1,'ConfidenceRegionDisplaySpacing',5)
opt=getoptions(h1);
opt.Title.String
opt.Title.String='';
setoptions(h1,opt)
figure(3),
sys=feedback(G_est,1);     %% get regions somehow.
h2=bodeplot(G_est);
setoptions(h2,'ConfidenceRegionNumberSD',5)
showConfidence(h2)
