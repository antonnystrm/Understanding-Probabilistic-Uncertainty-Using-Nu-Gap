close all
clear
clc
s=tf('s');
G_nom=2/(1+s/10)^3;             % 2->8 for instability in the mean

h1 = 0.1;
Gd = c2d(G_nom,h1);
frequencies = [0 0.05];
N=100;
u = idinput(N,'rbs',frequencies);
[y,t] = lsim(Gd,u);
w=linspace(0,75,1000);

n=1000;

dists=zeros(length(w),n);
models=[G_nom];
for k=1:n
    noiselevel = 0.1*std(y);
    e = noiselevel*randn(N,1);
    yn = y + e;

    noisydata = iddata(yn,u,h1);
    
    nz=0;np=3;                          %%model structure dependent
    G_est=tfest(noisydata,np, nz);
    dists(:,k)=chordal_dist(G_nom,G_est,w);
    models(k)=G_est;
    k
end
%%
figure(1)
surf(1:n,w,dists,'FaceAlpha',0.5, 'EdgeColor','none')
xlabel('Realization no.')
ylabel('Frequency')
zlabel('Chordal Distance')
title('Realizations of chordal distance to mean')
colorbar
%%
[M,I] = max(dists)
figure(2)
histogram(M);
titl='realizations of the \nu-gap';
title(sprintf('Histogram of %d %s \nbetween nominal model and realization.',n,titl));
xlabel('\nu-gap')
ylabel('Quantity')

figure(3)
histogram(w(I));
title(sprintf('Histogram of %d of peak chordal distance frequency.',n));
xlabel('Frequency')
ylabel('Quantity')

figure(4)
imagesc(1:n,w,dists)
colorbar
caxis([0 0.3]);
title(sprintf('Color map of %d realizations of point-wise \nchordal distance to the nominal model.',n));
xlabel('Realization no.')
ylabel('Frequency ')
%% plot the first m realizations on top of eachother & nyquist in same plot.
m=min(n,15);
figure(3)
hold on
figure(4)
hold on
for k=1:m
    figure(3)
    plot(w,dists(:,k),'b.')
    figure(4)
    nyquistplot(models(k),'b')
    k
end
nyquistplot(G_nom,'r')
