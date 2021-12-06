%% DMK33GR0134 camera : read image 
clearvars;

%% add relative path to current path
addpath('..\shared functions folder');
addpath('..\common subfunctions');
addpath('..\read and write files');
%% open image as tiff
[Filename,Foldername] = uigetfile('\\bazar\acousto-optique\Experimental Data\*.tiff','MultiSelect','on');

if iscell(Filename)==0
   Filename = {Filename};
end
Nfiles = length(Filename);

%% define a camera
MyDMK = camera('DMK33GR0134')    ;
MyDMK.format = 8;
MyDMK.wavelength = 1064e-9;
MyDMK.IntegrationTime = 60e-6;

%% import image
IM = double( importdata([Foldername,Filename{1}]) );
BG = double( importdata([Foldername,Filename{2}]) );

%% exemple 
% IM = double( importdata('D:\Data\Mai\2020-03-06\MAIN_J0_10us_14596.tiff') );

%% plot image

figure(1)
imagesc( 1e3*MyDMK.x_cam*MyDMK.dpixel , 1e3*MyDMK.y_cam*MyDMK.dpixel , IM - BG)
axis equal
cb = colorbar ;
%axis([1.7 2.6 0.8 2.9])
xlabel('mm')
ylabel('mm')
ylabel(cb,'Intensity in counts')
title('no pump')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% save as png

%% option code to retreive FWHM of gaussian function:

% find maximum position in the image:
[c,ind] = max(IM(:) - BG(:));
[row,col] = ind2sub(size(IM),ind);

fwhm_x = FWHM(IM(row,:) - BG(row,:),MyDMK.x_cam*MyDMK.dpixel)*1e3;
fwhm_y = FWHM(IM(:,col) - BG(:,col),MyDMK.y_cam*MyDMK.dpixel)*1e3;

Beam_x = 2/sqrt(-2*log(0.5))*fwhm_x
Beam_y = 2/sqrt(-2*log(0.5))*fwhm_y

figure;
plot(1e3*(MyDMK.x_cam*MyDMK.dpixel-MyDMK.x_cam(col)*MyDMK.dpixel),...
      IM(row,:) - BG(row,:));
hold on
plot(1e3*(MyDMK.y_cam*MyDMK.dpixel-MyDMK.x_cam(row)*MyDMK.dpixel),...
      IM(:,col) - BG(:,col));
legend('fwhm_x = ','fhwm_y = ')
%%


