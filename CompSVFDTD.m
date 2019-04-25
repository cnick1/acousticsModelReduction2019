% Comparing svd of the FDTD
close all
clear all
load SVFDTD

semilogx(sn4, 'bo')
hold on
semilogx(sn10, 'ro')
semilogx(sn20, 'go')
semilogx(sn40, 'ko')

xlabel('3n^2')

ylabel('Singular Values')
legend('n=4','n=10','n=20','n=40')

figure

plot(sn40, 'ko')
hold on
plot(sn40BC, 'ro')

xlabel('3n^2')
ylabel('Singular Values')
legend('No Buildings','One Building')