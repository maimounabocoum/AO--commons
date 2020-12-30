%% this VI reads .dat file recorded using labview VI
clearvars;
MU1 = importdata('D:\Data\Mai\2020-12-30\phase\MU1.dat');
MU2 = importdata('D:\Data\Mai\2020-12-30\phase\MU2.dat');

figure(1);
subplot(121)
imagesc(MU1)
title('\mu')
subplot(122)
imagesc(sqrt(MU2-MU1.^2))
title('\sigma')