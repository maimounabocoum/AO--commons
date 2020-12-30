%% add relative path to current path

addpath('..\shared functions folder');

%% load data and initialize classes
clearvars

%% load sequence
[Filename,Foldername] = uigetfile('D:\Data\Mai\2020-10-22\*.tif','MultiSelect','on');
if iscell(Filename)==0
   Filename = {Filename};
end
Nfiles = length(Filename);

%% load images of interest
% REF = double( importdata('Q:\datas\2019-12-10\EXP100\Ref_OnlyP40uW_nofiler.tif') );

BG = double( importdata([Foldername,Filename{2}]) );
IM = double( importdata([Foldername,Filename{1}]) );

% BG = double( importdata('D:\Data\Mai\2020-10-08\BG_input950mWW_Exposure50000us_16bit_100Hz_verrin0mmPlus17cm_contact.tif') );
% IM = double( importdata('D:\Data\Mai\2020-10-08\IM_input950mWW_Exposure50000us_16bit_100Hz_verrin0mmPlus17cm_contact.tif') );


% define a camera
MyXimea = camera('xiB-64');
MyXimea.format = 8;
MyXimea.wavelength = 780e-9;
MyXimea.IntegrationTime = 500e-6;
%MyXimea = MyXimea.ResizePixels(1024,1024);
% calibration (1/3.6604e-05)   : 29-09-2020
% calibration (4.9/1.0179e-04) :2010-10-02
% calibration (5.3/2.9673e-04) :2010-10-08
% get intensity and total power
% [Frame,P_tot]           = GetIntensity( IM - mean(mean(IM(1:10,1:10)))  , 0*BG , MyXimea );
% P_tot = trapz(MyXimea.x_cam*MyXimea.dpixel,trapz(MyXimea.y_cam*MyXimea.dpixel,(IM-BG) ));


NormFactor = ((1e-6)/MyXimea.IntegrationTime)*(5.3/2.9673e-04) ;
P_tot = trapz(MyXimea.x_cam*MyXimea.dpixel,trapz(MyXimea.y_cam*MyXimea.dpixel, NormFactor*(IM-BG) ));
mean(1e-4*NormFactor*(IM(:)-BG(:)))

figure(7)
imagesc( 1e3*MyXimea.x_cam*MyXimea.dpixel , 1e3*MyXimea.y_cam*MyXimea.dpixel , 1e-4*NormFactor*(IM-BG) )
cb = colorbar ;
caxis([0 0.1])
xlabel('mm')
ylabel('mm')
ylabel(cb,'Intensity in \mu W /cm^2')
title(['Total Power is ',num2str(P_tot),'\mu W - exposure 500\mu s'])


%%
pos =48-[48 40 30 20 10 0]+ 15 ;
fwhm = 2*[7.8 11.7 18 20 23 31];
figure(1)
subplot(121)
hold off
plot(pos,2*fwhm,'linewidth',4)
hold on
plot(pos,(0.78)./sin(atan(0.5./pos)),'linewidth',4,'marker','o')
xlabel('distance from iris (mm)')
ylabel('speckle fwhm (\mu)')
title('Evolution of speckle size (through 1mm iris)')
subplot(122)
hold off
plot(pos,[0.3327 0.1721 0.0942 0.0615 0.0456 0.0353]/0.3327,'linewidth',4,'marker','o')
hold on
NA = sin( atan( 0.5./(pos-15)) ) ;
plot( pos , NA/NA(1) , 'linewidth',4)
hold on
plot( pos , NA.^2/NA(1)^2,'linewidth',4)
xlabel('distance from iris (mm)')
ylabel('intensity (\mu W/cm^2)')
title('Evolution of collected Power')


%%
figure(7)
% subplot(221)
% imagesc( 1e3*MyXimea.x_cam*MyXimea.dpixel , 1e3*MyXimea.y_cam*MyXimea.dpixel , 1e-4*NormFactor*(IM-BG) )
% cb = colorbar ;
% xlabel('mm')
% ylabel('mm')
% ylabel(cb,'Intensity in \mu W /cm^2')
% title(['Total Power is ',num2str(P_tot),'\mu W - exposure 50ms'])
% subplot(222)
imagesc( 1e6*MyXimea.x_cam*MyXimea.dpixel-min(1e6*MyXimea.x_cam*MyXimea.dpixel ) , 1e6*MyXimea.y_cam*MyXimea.dpixel - min(1e6*MyXimea.y_cam*MyXimea.dpixel), 1e-4*NormFactor*(IM-BG) )
cb = colorbar ;
% caxis([0 0.2])
xlabel('\mu m')
ylabel('\mu m')
ylabel(cb,'Intensity in \mu W /cm^2')
title('Zoom')
% subplot(223)
% hist(1e-4*NormFactor*(IM(:)-BG(:)),100)
% xlabel('Intensity(\mu W /cm^2)')
% ylabel('number of pixels')
% title('hitogram of pixels')
% subplot(224)
% im = NormFactor*(IM(2500,:)-BG(2500,:)); 
% [ACorrIm,lags] = xcorr(im);
% imagesc(ACorrIm)
% plot(1e6*lags*MyXimea.dpixel ,ACorrIm,'o-');
% xlim([-250 250])
% xlabel('\mu m')



%%


















%%



