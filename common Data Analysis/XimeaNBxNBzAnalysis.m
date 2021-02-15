%%   this is a loop analysis routine %%%
%% add relative path to current path

addpath('..\shared functions folder');
addpath('..\common subfunctions')
addpath('..\read and write files');
addpath('D:\AO--simulations\Scan routines')

%% load data and initialize classes
clearvars -except P_tot_US  P_tot_OFF

%% import LogFile
[NB,nuX0,nuZ0] = ReadLogFile('D:\Data\Louis\2021-01-22\gel_retourne\sans_PJ\1-10'); % loads LogFile


%% load sequence corresponding to log file
[Filename,Foldername] = uigetfile('D:\Data\Louis\2021-01-22\gel_retourne\sans_PJ\1-10\*.tiff','MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end

%% load images of interest, BG and REF

IM = double( importdata([Foldername,Filename{10}]) );
BG      = double( importdata('D:\Data\Louis\2021-01-22\refs\BG.tiff') );
%MAIN    = double( importdata('D:\Data\Mai\2020-02-05\OBJ\OBJ_45268.tiff') );
%REF     = double( importdata('D:\Data\Mai\2020-02-05\REF\REF_48230.tiff') );

% figure;imagesc(BG);title('background');cb = colorbar;ylabel(cb,'counts');
% figure;imagesc(MAIN);title('main only');cb = colorbar;ylabel(cb,'counts');
% figure;imagesc(REF);title('ref only');cb = colorbar;ylabel(cb,'counts');

%% define a camera
MyXimea = camera('xiB-64')    ;
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 100e-6;
MyXimea = MyXimea.ResizePixels(size(BG,2),size(BG,1));

%% set references need for deconvolution
[Frame,~]               = GetIntensity( IM , BG , MyXimea );
% [Frame_ref,P_tot_ref]   = GetIntensity( REF , BG , MyXimea );
% [Frame_main,P_tot_main] = GetIntensity( MAIN , BG , MyXimea );

% figure; imagesc(100*Frame_ref);
% title(['ref P_{tot}(\mu W)=',num2str(1e6*P_tot_ref)]);cb = colorbar;ylabel(cb,'\mu W /cm^2');
% figure; imagesc(100*Frame_main);
% title(['main P_{tot}(\mu W)=',num2str(1e6*P_tot_main)]);cb = colorbar;ylabel(cb,'\mu W /cm^2');


%% Fourier class for data analysis

F = TF2D(MyXimea.Nx_cam , MyXimea.Ny_cam , 1/(MyXimea.dpixel) , 1/(MyXimea.dpixel) );
FrameFFT = F.fourier( Frame );
figure(1)
imagesc(abs(FrameFFT))
colorbar
caxis([0 0.5e-8])
myFilter = ImageFilter( [880 320 250 250] );
%myFilter = ImageFilter( [470 560 470 560] ); % center
BW = myFilter.getROI(F.Nx,F.Nz);
myFilter.DrawROI;
%%              -------------- begin loop -----------------%%

% Nfft = 2^11;
% G = TF2D( Nfft , Nfft , (Nfft-1)*nuX0 , (Nfft-1)*nuZ0 );
% ObjectFFT = zeros(Nfft , Nfft);



%for loop1 = 56%:length(NBloop)
%     loop_scan = find( NB(:,1)==NBloop(loop1) ) ;
%     for loop = loop_scan'%1:Nfiles
%     IndexRecord = ExtractIndex(Filename{loop}) - IndexRecord0 ;
%     Nbx = NB(mod(IndexRecord,Nfiles)+1,2) ;
%     Nbz = NB(mod(IndexRecord,Nfiles)+1,3) ;
%   ObjectFFT2((2^4+1)+Nbz,(2^4+1)+Nbx) = mean(mean(Gag.raw(5000:6000,loop_scan),1),2);
%     end
%   
% end
% figure;surfc(ObjectFFT2(18:27,12:22))

NBloop = unique(NB(:,1)); 
IndexRecord0    = ExtractIndex(Filename{1},'ximea'); % index of first frame
IM_rec          = zeros(F.Nz,F.Nx);
ImageCorr       = zeros(1,Nfiles);
P_tot           = zeros(1,Nfiles);

% loop_scan = find( NB(:,1)==NBloop(loop1-1) ) ;    
% loop_scan'
for loop = 1:Nfiles
    
%     if loop == 1
%     IM = double( importdata([Foldername,Filename{1}]) );
%     [Frame0,P_frame0] = GetIntensity( IM , BG , MyXimea ) ;     
%     end
    
    IM = double( importdata([Foldername,Filename{loop}]) );
    IndexRecord = ExtractIndex(Filename{loop},'ximea') - IndexRecord0 ; % Ximea program
    %IndexRecord = ExtractIndex(Filename{loop},'labview') - IndexRecord0 ;% Labview
    
    % get current image
    [Frame,P_frame] = GetIntensity( IM , BG , MyXimea );

    % filter out tagged photons
    FrameFFT      = F.fourier( Frame );
    % figure;imagesc(log(abs(FrameFFT))); myFilter.DrawROI;
    FilteredFrame = F.ifourier( FrameFFT.*BW) ;
       
    FilteredFrame = abs(FilteredFrame).^2;%./(Frame_ref);
    % figure;imagesc(FilteredFrame);
    
    %FilteredFrame(abs(FilteredFrame) > 10 ) = 0 ;
    %FilteredFrame(800:end,:)=0;
    % evaluation of total power
    temp = FilteredFrame;
    P_tot(loop) = sum( temp(:)*MyXimea.dpixel*MyXimea.dpixel );

    % indirect reconstruction
    
    Nbx     = NB(loop,2) ;
    Nbz     = NB(loop,3) ;
    PHASE   = NB(loop,4);
    
%     DecalZ  =   0.3;
%     s = exp(2i*pi*DecalZ*Nbz);
%     ObjectFFT((2^4+1)+Nbz,(2^4+1)+Nbx) = ObjectFFT((2^4+1)+Nbz,(2^4+1)+Nbx) + s*exp(1i*2*pi*PHASE)*P_tot(loop);
%     ObjectFFT((2^4+1)-Nbz,(2^4+1)-Nbx) = conj( ObjectFFT((2^4+1)+Nbz,(2^4+1)+Nbx));

    % direct reconstruction
    IM_rec = IM_rec + exp(1i*2*pi*PHASE)*FilteredFrame;
     
% figure(1)
% subplot(121)  
% imagesc( F.x*1e3 , F.z*1e3 , 100*FilteredFrame )
% cb = colorbar ;
% xlabel('mm')
% ylabel('mm')
% ylabel(cb,'Intensity in \mu W /cm^2')
% title(['P_{tot}',num2str(1e6*P_tot(loop)),'\mu W , (NBx,NBz,phase)=(',num2str(Nbx),',',num2str(Nbz),',',num2str(PHASE),')'])
% % ImageCorr(loop) = corr2(Frame,Frame0)    ;
% subplot(122)    
% % FilteredFrame_fft = F.fourier(FilteredFrame);
% % FilteredFrame_fft(512:514,512:514)= 0;
% imagesc( F.fx*1e3 , F.fz*1e3 , real( IM_rec ) )
% colorbar
% axis([-1 1 -1 1]*0.5e7)
% drawnow
%saveas(gcf,['D:\Data\Mai\2020-01-31\data analysis\PJ\',Filename{loop}],'png');
end

% figure(1)
% subplot(121)
% imagesc( F.x*1e3 , F.z*1e3 , real( IM_rec ) )
% caxis([-0.08 0.08])
% colorbar
% xlabel('x(mm)')
% ylabel('z(mm)')
% title('Real(4-phase Hologram)')
% subplot(122)
% imagesc( F.x*1e3 , F.z*1e3 , imag( IM_rec ) )
% caxis([-0.08 0.08])
% colorbar
% title('Imag(4-phase Hologram)')
% xlabel('x(mm)')
% ylabel('z(mm)')
% 
% saveas(gcf,['D:\Data\Mai\2020-01-31\data analysis\PJ\image',num2str(loop1)],'png');

% cal = 2.75;
% nuX = Nbx*nuX0;
% nuZ = Nbz*nuZ0;
% alpha = nuX/nuZ;
% 
% calX = 6.6;
% calY = 7;
% 
% [X,Z] = meshgrid(F.x,F.z);
% PHASE0 = exp(2*1i*pi*(calX*nuX*X +calY*nuZ*Z));
% PHASE0 = PHASE0(100:750,100:800);
% 
% figure(2)
% subplot(121); imagesc(real(IM_rec(100:750,100:800))); colorbar ;
% title(['Real : (NBx,NBz)=(',num2str(Nbx),',',num2str(Nbz),')'])
% subplot(122); imagesc(imag(IM_rec(100:750,100:800))); colorbar ;
% title(['Imag : (NBx,NBz)=(',num2str(Nbx),',',num2str(Nbz),')'])
%subplot(223); imagesc(abs(IM_rec(100:750,100:800))); colorbar ;
% subplot(223); imagesc(cal*F.x*1e3,cal*F.z*1e3,angle(PHASE0)); colorbar ;
% %subplot(224); imagesc(cal*F.x*1e3,cal*F.z*1e3,angle(IM_rec(100:750,100:800).*PHASE0)); colorbar ;
% subplot(224); imagesc(cal*F.x*1e3,cal*F.z*1e3,abs(IM_rec(100:750,100:800))); colorbar ;

% r = IM_rec(100:750,100:800).*PHASE0;
% r = r(200:500,200:500);

% ObjectFFT((2^4+1)+Nbz,(2^4+1)+Nbx) = (mean(abs(r(:))));

% mean(Gag.raw(5000:6000,:),2)

%caxis([0 7e-3])
% title(['Real : (NBx,NBz)=(',num2str(Nbx),',',num2str(Nbz),')'])
% saveas(gcf,['D:\Data\Mai\2020-01-14\data analysis\Imag(NBx,NBz)(',num2str(Nbx),',',num2str(Nbz),')'],'png');

% end

%%

% [Xf,Zf] = meshgrid(-5:5,1:10);
% 
% PhasL =  1.1239*Xf + 0.718*Zf - 3.939;
% figure(6);
% subplot(121); surfc( -PhasL + UnwrapPHASE(ObjectFFT(18:27,12:22),5,5) );
% xlabel('NbX')
% ylabel('NbZ')
% view([-58.8340 17.8721])
% colorbar
% subplot(122); imagesc(PhasL);
% colorbar
% caxis([])
% xlabel('NBx')
% ylabel('NbZ')

%% interpolation of main object onto this grid:
% MAIN_FFT = F.fourier(MAIN);
% cal = 7.5;
% MAIN_FFT( F.Nz/2+1,:) = 0;
% MAIN_FFT( (F.Nz/2+1)+max(NB(:,3)):end,:) = 0;
% MAIN_FFT( 1:((F.Nz/2+1)-max(NB(:,3))),:) = 0;
% MAIN_FFT(:,(F.Nx/2+1)+max(NB(:,2)):end) = 0;
% MAIN_FFT(:,1:((F.Nx/2+1)-max(NB(:,2)))) = 0;
% MAIN_FFT(:,1:((F.Nx/2+1)-max(NB(:,2)))) = 0;
% %MAIN_FFT([F.Nz/2,F.Nz/2+2],F.Nx/2+1) = MAIN_FFT([F.Nz/2,F.Nz/2+2],F.Nx/2+1);
% figure(4);imagesc((1/cal)*F.fx/G.dfx,(1/cal)*F.fz/G.dfz,abs(MAIN_FFT))
% axis([-15 15 -10 10])
% xlabel('NbX')
% ylabel('NbZ')
% colorbar
% 
% MAIN_filtered = F.ifourier(MAIN_FFT);
% MAIN_filtered = MAIN_filtered - ones(Nfft,1)*MAIN_filtered(1,:);
% figure(5);imagesc(cal*F.x*1e3,cal*F.z*1e3,MAIN_filtered)

%% loop to get fourier componant:
% load .dat
% DATA = importdata('D:\Data\Louis\2021-01-07\premier image\avecLens.dat');
% P_tot = DATA(:,1);
% 
%  NB(1:4,:) = [];
%  P_tot(1:4) = [];

Nfft = 2^10;
G0 = TF2D( Nfft , Nfft , (Nfft-1)*nuX0 , (Nfft-1)*nuZ0 );
ObjectFFT = zeros(Nfft , Nfft);

for loop = 5:length(P_tot)
    Nbx     = NB(loop,2) ;
    Nbz     = NB(loop,3) ;
    PHASE   = NB(loop,4);
    
    DecalZ  =   0; % ??
    DecalX  =   0; % ??
    s =  1i*exp(2i*pi*(DecalZ*Nbz + DecalX*Nbx));
   
    ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) = ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) + s*exp(1i*2*pi*PHASE)*P_tot(loop);
    ObjectFFT((Nfft/2+1)-Nbz,(Nfft/2+1)-Nbx) = conj( ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) );   

end

%% inversion using simulated Matrix
[U,S,V] = svd(G.M);
 Splus = 0*S;
 seuil = S(200,200);
 Splus(S>seuil) = 1./S(S>seuil);
 Splus = Splus';
 Minv = V*Splus*U';
 Ireconst = Minv*P_tot(5:end)';

 Ireconst = squeeze(reshape( Ireconst ,[1,G.Nx,G.Nz]))';
 figure;imagesc(G.x,G.z,real( Ireconst));
 colorbar;
% axis([-15 15 -10 10])
% %caxis([-2 2]*10e-8)
% xlabel('NbX')
% ylabel('NbZ')

%%

% get FFT

Reconstruct = G0.ifourier( ObjectFFT );
% % I = ifft2(ifftshift(ObjectFFT));
% Reconstruct = Reconstruct - ones(Nfft,1)*Reconstruct(1,:);
% % I = ifftshift(I,2);
figure(2);
imagesc(G0.x*1e3,G0.z*1e3,real(Reconstruct))
title('reconstructed AO image')
xlabel('x(mm)')
ylabel('z(mm)')
cb = colorbar;
ylabel(cb,'a.u.')

% l = 3.1746/2;
% rectangle('Position',[-1+2*l,-4.8,l,5*l],...
%           'Curvature',[0,0],...
%          'LineWidth',2)
% rectangle('Position',[-1,-4.8,l,5*l],...
%           'Curvature',[0,0],...
%          'LineWidth',2)
% rectangle('Position',[-1-2*l,-4.8,l,5*l],...
%           'Curvature',[0,0],...
%          'LineWidth',2)


    