%%   this is a loop analysis routine %%%
%% add relative path to current path

addpath('..\shared functions folder');
addpath('..\common subfunctions')
addpath('..\read and write files');

%% load data and initialize classes
clearvars

%% import LogFile
[NB,nuX0,nuZ0] = ReadLogFile('');
NB(1:4,:) = [];

%% load sequence
[Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end

%% load images of interest

%REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_nofiler.tif') );
BG      = double( importdata('D:\Data\Mai\2020-01-13\NbXNbZ\GB\BG_48524.tiff') );
MAIN    = double( importdata('D:\Data\Mai\2020-01-14\MAIN\MAIN_15181.tiff') );
REF     = double( importdata('D:\Data\Mai\2020-01-13\REF\ref_12911.tiff') );

% figure;imagesc(BG);title('background');cb = colorbar;ylabel(cb,'counts');
% figure;imagesc(MAIN);title('main only');cb = colorbar;ylabel(cb,'counts');
% figure;imagesc(REF);title('ref only');cb = colorbar;ylabel(cb,'counts');

%% define a camera
MyXimea = camera('xiB-64')    ;
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 100e-6;
MyXimea = MyXimea.ResizePixels(1024,1024);

%% set references need for deconvolution
[Frame_ref,P_tot_ref]   = GetIntensity( REF , BG , MyXimea );
[Frame_main,P_tot_main] = GetIntensity( MAIN , BG , MyXimea );

figure;imagesc(100*Frame_ref);title(['ref P_{tot}(\mu W)=',num2str(1e6*P_tot_ref)]);cb = colorbar;ylabel(cb,'\mu W /cm^2');
figure;imagesc(100*Frame_main);title(['main P_{tot}(\mu W)=',num2str(1e6*P_tot_main)]);cb = colorbar;ylabel(cb,'\mu W /cm^2');


%% fourier class
N   = 2^10;
F = TF2D( N , N , 1/(MyXimea.dpixel) , 1/(MyXimea.dpixel) );
%% define filter for FFT in pixel
myFilter = ImageFilter( [600 320 90 90] );
%myFilter = ImageFilter( [470 560 470 560] ); % center
BW = myFilter.getROI(N,N);
myFilter.DrawROI;
%% -------------- begin loop 
IndexRecord0 = ExtractIndex(Filename{1});
IM_rec  = zeros(N,N);

ImageCorr = zeros(1,Nfiles);
G = TF2D( 2^5 , 2^5 , (2^5-1)*nuX0 , (2^5-1)*nuZ0 );
ObjectFFT = zeros(2^5 , 2^5);


for loop = 1:Nfiles
    
%     if loop == 1
%     IM = double( importdata([Foldername,Filename{1}]) );
%     [Frame0,P_frame0] = GetIntensity( IM , BG , MyXimea ) ;     
%     end
    
    IM = double( importdata([Foldername,Filename{loop}]) );
    IndexRecord = ExtractIndex(Filename{loop}) - IndexRecord0 ;
    [Frame,P_frame] = GetIntensity( IM , BG , MyXimea );
    
    
    % filter out tagged photons
    FrameFFT      = F.fourier( Frame );
    % figure;imagesc(log(abs(FrameFFT)));
    FilteredFrame = F.ifourier( FrameFFT.*BW) ;
    FilteredFrame = 2*abs(FilteredFrame).^2;%./(Frame_ref);
    %FilteredFrame(abs(FilteredFrame) > 10 ) = 0 ;
    %FilteredFrame(800:end,:)=0;
    % evaluation of total power
    temp = FilteredFrame;
    P_tot(loop) = sum( temp(:)*MyXimea.dpixel*MyXimea.dpixel );
    

    
    % indirect reconstruction
    
    Nbx = NB(mod(IndexRecord,Nfiles)+1,2) ;
    Nbz = NB(mod(IndexRecord,Nfiles)+1,3) ;
    PHASE= NB(mod(IndexRecord,Nfiles)+1,4)+pi/2;
    
    DecalZ  =   0.3;
    s = exp(2i*pi*DecalZ*Nbz);
    ObjectFFT((2^4+1)+Nbz,(2^4+1)+Nbx) = ObjectFFT((2^4+1)+Nbz,(2^4+1)+Nbx) + s*exp(1i*2*pi*PHASE)*P_tot(loop);
    ObjectFFT((2^4+1)-Nbz,(2^4+1)-Nbx) = conj( ObjectFFT((2^4+1)+Nbz,(2^4+1)+Nbx));

    % direct reconstruction
    IM_rec = IM_rec + exp(1i*2*pi*PHASE)*FilteredFrame;
%     
%     figure
% imagesc( F.x*1e3 , F.z*1e3 , 100*FilteredFrame )
% cb = colorbar ;
% xlabel('mm')
% ylabel('mm')
% ylabel(cb,'Intensity in \mu W /cm^2')
% title(['P_{tot}',num2str(1e6*P_tot(loop)),'\mu W , (NBx,NBz,phase)=(',num2str(Nbx),',',num2str(Nbz),',',num2str(PHASE),')'])
% % ImageCorr(loop) = corr2(Frame,Frame0)    ;
% % subplot(122)    
% % FilteredFrame_fft = F.fourier(FilteredFrame);
% % %FilteredFrame_fft(512:514,512:514)= 0;
% % imagesc( F.fx*1e3 , F.fz*1e3 , UnwrapPHASE(angle(FilteredFrame_fft),512,512))
% % colorbar
% % axis([-1 1 -1 1]*0.5e7)
% % caxis([-30 30])
% drawnow
% saveas(gcf,['D:\Data\Mai\2020-01-13\4-phase\',Filename{loop}],'png');
end


% figure(2)
% subplot(131); imagesc(real(IM_rec(100:750,100:800))); colorbar ;
% title(['Real : (NBx,NBz)=(',num2str(Nbx),',',num2str(Nbz),')'])
% subplot(132); imagesc(imag(IM_rec(100:750,100:800))); colorbar ;
% title(['Imag : (NBx,NBz)=(',num2str(Nbx),',',num2str(Nbz),')'])
% subplot(133); imagesc(abs(IM_rec(100:750,100:800))); colorbar ;
% caxis([0 7e-3])
% title(['Real : (NBx,NBz)=(',num2str(Nbx),',',num2str(Nbz),')'])
% saveas(gcf,['D:\Data\Mai\2020-01-14\data analysis\Imag(NBx,NBz)(',num2str(Nbx),',',num2str(Nbz),')'],'png');

%% interpolation of main object onto this grid:
MAIN_FFT = F.fourier(MAIN);
cal = 7.5;
MAIN_FFT(F.Nz/2+1,:) = 0;
MAIN_FFT((F.Nz/2+1)+6:end,:) = 0;
MAIN_FFT(1:((F.Nz/2+1)-6),:) = 0;
MAIN_FFT(:,(F.Nx/2+1)+6:end) = 0;
MAIN_FFT(:,1:((F.Nx/2+1)-6)) = 0;
MAIN_FFT(:,1:((F.Nx/2+1)-6)) = 0;
MAIN_FFT([F.Nz/2,F.Nz/2+2],F.Nx/2+1) = 0.01*MAIN_FFT([F.Nz/2,F.Nz/2+2],F.Nx/2+1);
figure(4);imagesc((1/cal)*F.fx/G.dfx,(1/cal)*F.fz/G.dfz,real(MAIN_FFT))
axis([-15 15 -10 10])
xlabel('NbX')
ylabel('NbZ')
colorbar
MAIN_filtered = F.ifourier(MAIN_FFT);
figure(5);imagesc(cal*F.x*1e3,cal*F.z*1e3,MAIN_filtered)

%% loop to get fourier componant:
Nfft = 2^10;

G = TF2D( Nfft , Nfft , (Nfft-1)*nuX0 , (Nfft-1)*nuZ0 );
ObjectFFT = zeros(Nfft , Nfft);

for loop = 1:length(P_tot)
    Nbx     = NB(loop,2) ;
    Nbz     = NB(loop,3) ;
    PHASE   = NB(loop,4);
    
    DecalZ  =   -0.2; % ??
    chirpZ  =   0.3;
    DecalX  =   -0.2; % ??
    s = exp(2i*pi*(DecalZ*Nbz + chirpZ*Nbz.^2 + DecalX*Nbx));
   
    ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) = ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) + s*exp(1i*2*pi*PHASE)*P_tot(loop);
    ObjectFFT((Nfft/2+1)-Nbz,(Nfft/2+1)-Nbx) = conj( ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) );   

end

% figure;imagesc(G.fx/G.dfx,G.fz/G.dfz,abs(ObjectFFT)), colorbar;
% axis([-15 15 -10 10])
% %caxis([-2 2]*10e-8)
% xlabel('NbX')
% ylabel('NbZ')



% get FFT

Reconstruct = G.ifourier( ObjectFFT );
% % I = ifft2(ifftshift(ObjectFFT));
% % I = I - ones(2^5,1)*I(1,:);
% % I = ifftshift(I,2);
figure(1);imagesc(G.x*1e3,G.z*1e3,real(Reconstruct))
title('reconstructed object')
xlabel('mm')
ylabel('mm')
cb = colorbar;


    