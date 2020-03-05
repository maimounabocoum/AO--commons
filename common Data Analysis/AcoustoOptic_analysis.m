%% common analysis for Phorefractive experiement 
clearvars;
close all ;

addpath('..\shared functions folder')

%% data loaded
% 2017-10-17 : JIONC presentation 

%% data analysis

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
   pitch = 0.2;
   [I,X,Z] = Reconstruct(NbX , NbZ, ...
                         NUX , NUZ ,...
                         x , z , ...
                         Datas - mean(mean(Datas(1:10,1:10))) , ...
                         SampleRate , 20 , c , pitch);      
        
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












