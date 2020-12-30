%%   this is a loop analysis routine %%%

%% add relative path to current path

addpath('..\shared functions folder');

%% load data and initialize classes
clearvars

%% load sequence
[Filename,Foldername] = uigetfile('D:\Data\Mai\2020-11-02*.tiff','MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end

%% load images of interest

%REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_nofiler.tif') );
BG = double( importdata('D:\Data\Mai\2020-11-02\references\BG_exp100us.tiff') );
MAIN = double( importdata('D:\Data\Mai\2020-11-02\references\IM1_100us.tiff') );
% REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_EXP100us.tif') );

%% define a camera
MyXimea = camera('xiB-64')    ;
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 100e-6;
MyXimea = MyXimea.ResizePixels(1024,1024);

%% set references need for deconvolution
% [Frame_ref,P_tot_ref]   = GetIntensity( REF , BG , MyXimea );
% [Frame_main,P_tot_main] = GetIntensity( MAIN , BG , MyXimea );


%% fourier class
N   = 2^10;
F = TF2D( N , N , 1/(MyXimea.dpixel) , 1/(MyXimea.dpixel) );
%% define filter for FFT in pixel
myFilter1 = ImageFilter( [653 531 95 95] ); % [centerX centerY Lx Ly]
myFilter2 = ImageFilter( [50 50 90 90] ); % center
BW1 = myFilter1.getROI(N,N);
BW2 = myFilter2.getROI(N,N);
% 
%  figure(1);imagesc(log(abs(FrameFFT)))Frame
% figure(1);imagesc(log(abs(Frame)))
%   myFilter1.DrawROI ;
%   myFilter2.DrawROI ;

%% -------------- begin loop 
% MU  = zeros(N,N);
% MU2 = zeros(N,N);
% P1 = zeros(1,Nfiles);
% P2 = zeros(1,Nfiles);
% 
% mu =  0;
% mu2 = 0;
% ImageCorr = zeros(1,Nfiles);

for loop = 1:Nfiles
    
%     if loop == 1
%     IM = double( importdata([Foldername,Filename{loop}]) );
%     [Frame0,P_frame0] = GetIntensity( IM , BG , MyXimea ) ;     
%     end
    
    IM = double( importdata([Foldername,Filename{loop}]) );
    [Frame,P_frame] = GetIntensity( IM , BG , MyXimea );
    
    % filter out tagged photons
%     FrameFFT = F.fourier( Frame );
%     FilteredFrame1 = F.ifourier( FrameFFT.*BW1) ;
%     FilteredFrame1 = 2*abs(FilteredFrame1).^2./(Frame_ref);
%     FilteredFrame2 = F.ifourier( FrameFFT.*BW2) ;
%     FilteredFrame2 = 2*abs(FilteredFrame2).^2./(Frame_ref);    
    
    % figure; imagesc(log(abs(FrameFFT)))
    % myFilter1.DrawROI ;
    % myFilter2.DrawROI ;
    % evaluation of total power

    
    Power(loop) = P_frame ;
 
    % ImageCorr(loop) = corr2(Frame,Frame0)    ;
end

%% plot result
    [Datas_mu1,Datas_std1,Datas_mu2,Datas_std2] = AverageDataBothWays( raw/(0.45*1e5) );
    t1 = (1/Fs2)*(1:length(Datas_mu1));
    t2 = (1/Fs1)*(1:length(Datas_mu2));
    
figure(2);
Iextract = 10:100;
plot(t1(Iextract)*1e3,1e6*Power(Iextract)-0.8);
hold on
plot(t1(Iextract)*1e3,1e6*Datas_mu1(Iextract))
%ylim([20.4 21])

CC = corrcoef( Power(Iextract), Datas_mu1(Iextract) ) ;
legend('CCD ximea','Photodiode')
title(['amplitude modulation 10mV - correlation ',num2str(CC(1,2))])
xlabel('time (ms)')
ylabel('Power (uW)')



%%






