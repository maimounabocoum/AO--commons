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
BG = double( importdata('D:\Data\Mai\2020-01-03\BG_51211.tiff') );
MAIN = double( importdata('D:\Data\Mai\2020-01-03\MAIN_51719.tiff') );
REF = double( importdata('D:\Data\Mai\2020-01-03\ref_51473.tiff') );

%% define a camera
MyXimea = camera('xiB-64')    ;
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 100e-6;
MyXimea = MyXimea.ResizePixels(1024,1024);

% get intensity and total power
[Frame,P_tot]           = GetIntensity( IM , BG , MyXimea );
[Frame_ref,P_tot_ref]   = GetIntensity( REF , BG , MyXimea );
[Frame_main,P_tot_main] = GetIntensity( MAIN , BG , MyXimea );

figure(1)
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
% myFilter = ImageFilter( [610 700 450 600] );
myFilter = ImageFilter( [410 720 120 120] );
BW       = myFilter.getROI(2^10,2^10);

figure(2)
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

%%
% inverse FFT
FilteredFrame = F.ifourier( FrameFFT.*BW) ;
FilteredFrame = 2*abs(FilteredFrame).^2;%./(Frame_ref) ;
% total sum in power
P_tot = sum( abs(FilteredFrame(:))*MyXimea.dpixel*MyXimea.dpixel );

figure(3)
imagesc( F.x*1e3 , F.z*1e3 , 100*FilteredFrame )
cb = colorbar ;
%caxis([15 115])
xlabel('mm')
ylabel('mm')
ylabel(cb,'Intensity in \mu W /cm^2')
title(['Total Power is ',num2str(1e6*P_tot),'\mu W'])
