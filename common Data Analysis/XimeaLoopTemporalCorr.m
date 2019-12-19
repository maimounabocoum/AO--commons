%%   this is a loop analysis routine %%%

%% add relative path to current path

addpath('..\shared functions folder');

%% load data and initialize classes
clearvars

%% load sequence
[Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end

%% load images of interest

%REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_nofiler.tif') );
BG = double( importdata('Q:\datas\2019-12-10\EXP100\BG_EXP100us.tif') );
MAIN = double( importdata('Q:\datas\2019-12-10\EXP100\Main_OnlyP9uW_EXP100us.tif') );
REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_EXP100us.tif') );

%% define a camera
MyXimea = camera('xiB-64')    ;
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 100e-6;

%% set references need for deconvolution
[Frame_ref,P_tot_ref] = GetIntensity( REF , BG , MyXimea );
[Frame_main,P_tot_main] = GetIntensity( MAIN , BG , MyXimea );


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
MU  = zeros(N,N);
MU2 = zeros(N,N);
P1 = zeros(1,Nfiles);
P2 = zeros(1,Nfiles);

mu =  0;
mu2 = 0;
ImageCorr = zeros(1,Nfiles);

for loop = 1:Nfiles
    
%     if loop == 1
%     IM = double( importdata([Foldername,Filename{loop}]) );
%     [Frame0,P_frame0] = GetIntensity( IM , BG , MyXimea ) ;     
%     end
    
    IM = double( importdata([Foldername,Filename{loop}]) );
    [Frame,P_frame] = GetIntensity( IM , BG , MyXimea );
    
    % filter out tagged photons
    FrameFFT = F.fourier( Frame );
    FilteredFrame1 = F.ifourier( FrameFFT.*BW1) ;
    FilteredFrame1 = 2*abs(FilteredFrame1).^2./(Frame_ref);
    FilteredFrame2 = F.ifourier( FrameFFT.*BW2) ;
    FilteredFrame2 = 2*abs(FilteredFrame2).^2./(Frame_ref);    
    
    % figure; imagesc(log(abs(FrameFFT)))
    % myFilter1.DrawROI ;
    % myFilter2.DrawROI ;
    % evaluation of total power
    frame1 = FilteredFrame1 ;
    frame2 = Frame.*BW2 ;
    P1(loop) = sum( abs(frame1(:))*MyXimea.dpixel*MyXimea.dpixel );
    P2(loop) = sum( abs(frame2(:))*MyXimea.dpixel*MyXimea.dpixel );
 
    % ImageCorr(loop) = corr2(Frame,Frame0)    ;
end

%% plot result
figure(2);plot(P1-mean(P1)); hold on ;plot(P2-mean(P2))
title( ['Corr = ', num2str(corr2(P1,P2)) ] )

