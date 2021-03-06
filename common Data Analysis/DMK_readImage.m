%% DMK33GR0134 camera : read image 
clearvars;

%% add relative path to current path
addpath('..\shared functions folder');
addpath('..\read and write files');
%% open image as tiff
[Filename,Foldername] = uigetfile('D:\Data\Mai\2020-09-02\*.tiff','MultiSelect','on');

if iscell(Filename)==0
   Filename = {Filename};
end
Nfiles = length(Filename);

%% define a camera
MyDMK = camera('DMK33GR0134')    ;
MyDMK.format = 8;
MyDMK.wavelength = 1064e-9;
MyDMK.IntegrationTime = 100e-6;

%% import image
IM = double( importdata([Foldername,Filename{2}]) );
BG = double( importdata([Foldername,Filename{1}]) );

%% exemple 
% IM = double( importdata('D:\Data\Mai\2020-03-06\MAIN_J0_10us_14596.tiff') );

%% plot image

figure(1)
%imagesc( 1e3*MyDMK.x_cam*MyDMK.dpixel , 1e3*MyDMK.y_cam*MyDMK.dpixel , IM - BG)
imagesc( MyDMK.x_cam*0.0047 , MyDMK.y_cam*0.0047 , IM - BG)
axis equal
cb = colorbar ;
%axis([1.7 2.6 0.8 2.9])
xlabel('mm')
ylabel('mm')
ylabel(cb,'Intensity in counts')
title('no pump')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% save as png




