%% add relative path to current path

addpath('..\shared functions folder');

%% load data and initialize classes
clearvars

%% load sequence
[Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
if iscell(Filename)==0
   Filename = {Filename};
end
Nfiles = length(Filename);
%% load images of interest

%REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_nofiler.tif') );
IM = double( importdata([Foldername,Filename{2}]) );
BG = double( importdata('D:\Data\Mai\2021-01-04\refs\BG.tiff') );
MAIN = double( importdata('D:\Data\Mai\2021-01-04\refs\MAIN.tiff') );
REF = double( importdata('D:\Data\Mai\2021-01-04\refs\REF.tiff') );

%% define a camera
MyXimea = camera('xiB-64')    ;
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 100e-6;
MyXimea = MyXimea.ResizePixels(1024,1024);

% get intensity and total power
[Frame,P_tot]           = GetIntensity( IM   , BG , MyXimea );
[Frame_ref,P_tot_ref]   = GetIntensity( REF  , BG , MyXimea );
[Frame_main,P_tot_main] = GetIntensity( MAIN , BG , MyXimea );

figure(11)
subplot(131)
imagesc( 1e3*MyXimea.x_cam*MyXimea.dpixel , 1e3*MyXimea.x_cam*MyXimea.dpixel , 100*Frame)
cb = colorbar ;
xlabel('mm')
ylabel('mm')
ylabel(cb,'Intensity in \mu W /cm^2')
title(['Total Power is ',num2str(1e6*P_tot),'\mu W'])
subplot(132)
imagesc( 1e3*MyXimea.x_cam*MyXimea.dpixel , 1e3*MyXimea.x_cam*MyXimea.dpixel , 100*Frame_ref)
cb = colorbar ;
xlabel('mm')
ylabel('mm')
ylabel(cb,'Intensity in \mu W /cm^2')
title(['Total Power is ',num2str(1e6*P_tot_ref),'\mu W'])
subplot(133)
imagesc( 1e3*MyXimea.x_cam*MyXimea.dpixel , 1e3*MyXimea.x_cam*MyXimea.dpixel , 100*Frame_main)
cb = colorbar ;
xlabel('mm')
ylabel('mm')
ylabel(cb,'Intensity in \mu W /cm^2')
title(['Total Power is ',num2str(1e6*P_tot_main),'\mu W'])



%% Definition of Fourier transform

F = TF2D( 2^10 , 2^10 , 1/(MyXimea.dpixel) , 1/(MyXimea.dpixel) );
FrameFFT = F.fourier( Frame );

% define filter in pixel
% myFilter = ImageFilter( [750 640 120 120] );
% myFilter = ImageFilter( [745 680 120 120] );
myFilter = ImageFilter( [790 535 165 115] );
BW       = myFilter.getROI(2^10,2^10);

figure(1)
%imagesc( 1e3*MyXimea.x_cam*MyXimea.dpixel , 1e3*MyXimea.x_cam*MyXimea.dpixel ,log(abs(FrameFFT)) )
imagesc(log(abs(FrameFFT)))
cb = colorbar ;
xlabel('mm')
ylabel('mm')
ylabel(cb,'Intensity in \mu W /cm^2')
title(['Total Power is ',num2str(1e6*P_tot),'\mu W'])
axis equal

% plot ROI (imagesc should be in pixels) in gcf
myFilter.DrawROI ;

%% loop on data 
[Frame_ref,P_tot_ref]   = GetIntensity( REF , BG , MyXimea );
[Frame_main,P_tot_main] = GetIntensity( MAIN , BG , MyXimea );

P_tot = [];
for loop = 1:Nfiles   
% inverse FFT
IM = double( importdata([Foldername,Filename{loop}]) );
[Frame,~]  = GetIntensity( IM , BG , MyXimea );
%FrameFFT = F.fourier( Frame );
%FilteredFrame = F.ifourier( FrameFFT.*BW ) ;
%FilteredFrame = abs(FilteredFrame).^2 ; %/(1.1+Frame_ref)

% total sum in power
P_tot(loop) = sum( abs(FilteredFrame(:))*MyXimea.dpixel*MyXimea.dpixel );

figure(2)
imagesc( F.x*1e3 , F.z*1e3 , FilteredFrame )
cb = colorbar ;
%caxis([0 10])
xlabel('mm')
ylabel('mm')
ylabel(cb,'Intensity in \mu W /cm^2')
title(['Total Power is ',num2str(1e6*P_tot(loop)),'\mu W'])
%
end
%% plot Ptot
    figure(3)
    hold off
    plot(repmat(2,1,length(P_tot)),100*P_tot/P_tot_main,'o')
    xlabel('Volt')
    ylabel('efficiency[%]')
    title('measured power')
%     
    %%

