%% this VI reads .dat file recorded using labview VI
clearvars;
MU1 = importdata('D:\Data\Mai\2021-01-04\DoublingEffect\MU1.dat');
MU2 = importdata('D:\Data\Mai\2021-01-04\DoublingEffect\MU2.dat');

figure(1);
subplot(121)
imagesc(MU1)
colorbar
title('\mu')
subplot(122)
imagesc(sqrt(MU2-MU1.^2))
colorbar
title('\sigma')

%% import LogFile
clearvars
[NB,nuX0,nuZ0] = ReadLogFile('D:\Data\Louis\2021-01-07\premier image\LogFile'); % loads LogFile
NB(1:4,:) = [];


DATA = importdata('D:\Data\Louis\2021-01-07\premier image\avecLens.dat');
P_tot = DATA(:,1);

for loop = 1:length(P_tot)
    Nbx     = NB(loop,2) ;
    Nbz     = NB(loop,3) ;
    PHASE   = NB(loop,4);
    
    DecalZ  =   -0.5; % ??
    DecalX  =   -0.2; % ??
    s =  1i*exp(2i*pi*(DecalZ*Nbz + DecalX*Nbx));
   
    ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) = ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) + s*exp(1i*2*pi*PHASE)*P_tot(loop);
    ObjectFFT((Nfft/2+1)-Nbz,(Nfft/2+1)-Nbx) = conj( ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) );   

end

figure(1);
imagesc(DATA)
colorbar
title('\mu')
%%