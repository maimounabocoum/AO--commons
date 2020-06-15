%% common analysis for Phorefractive experiement 
clearvars;
close all ;

addpath('..\shared functions folder')
load('Q:\datas\experimental datas\2017-10-17\AgarSL102_TypeOfSequence_JM_NbZ_8_NbZ_10_15h25_15.mat')
%% data loaded
% 2017-10-17 : JIONC presentation 

%% data analysis


% define main frequency ( older version of data )
DurationWaveform = 20;
Nbtot = 192;
pitch = 0.2;
nuZ0 = (1e6)/(c*DurationWaveform*1e3);       % Pas fréquence spatiale en Z (en mm-1)
nuX0 = 1/(Nbtot*pitch);                      % Pas fréquence spatiale en X (en mm-1)


% edit colormap
D=[1 1 1;0 0 1;0 1 0;1 1 0;1 0 0;];
F=[0 0.25 0.5 0.75 1];
G=linspace(0,1,256);
cmap=interp1(F,D,G);

switch TypeOfSequence
    
    case 'OF'
        
    Hfinal = figure(1) ;
    imagesc(x,z,Datas - mean(mean(Datas(1:10,1:10))) )
    xlabel('x(mm)')
    ylabel('z(mm)')
    cb = colorbar;
    colormap(parula)
    ylabel(cb,'a.u')
    title('Ondes Focalisées, \Delta t = 1.4 \mu s')
    set(findall(Hfinal,'-property','FontSize'),'FontSize',15)
       
    case 'JM'
    
        % inversion from FFT
   
   [I,X,Z] = Reconstruct(NbX , NbZ, ...
                         NUX , NUZ ,...
                         x , z , ...
                         Datas - mean(mean(Datas(1:10,1:10))) , ...
                         SampleRate , DurationWaveform , c , nuX0 , nuZ0);      
        
    Hfinal = figure(1) ;
    imagesc(x,z, I)
    xlabel('x(mm)')
    ylabel('z(mm)')
    cb = colorbar;
    colormap(parula)
    ylabel(cb,'a.u')
    title('Ondes Focalisées, \Delta t = 1.4 \mu s')
    set(findall(Hfinal,'-property','FontSize'),'FontSize',15)
        
        
end

%%












