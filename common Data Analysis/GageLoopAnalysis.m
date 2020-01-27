%% gage loop analysis
%clearvars
%load('Q:\datas\2019-12-10\Ref_100Hz_41uW_filter_main2uW_RepHz_100_16h27_53.mat');
load('D:\Data\Mai\2020-01-13\gage');
%% =============== add subfolders ===============  %%

%  addpath('Q:\AO--commons\Plottings')
%  addpath('..\..\AO--commons\shared functions folder')
%  addpath('gui');
%  addpath('sequences');
%  addpath('subfunctions');
 
 %% =============== load datas ===============  %%
 
 % data : 13/11/2019:
 %load('Q:\datas\2019-12-10\Ref_100Hz_41uW_filter_main2uW_RepHz_100_16h27_53.mat')
 
 
 %% =============== screen out datas ===============  %%

     h = 6.6e-34;
     lambda = 780e-9;
     Ephoton = h*(3e8/lambda);
     Fs1 = SampleRate*1e6 ; % gage sampling frequency
     Fs2 = 100 ;            % triggerbox sampling frequency
     % definition of 2 structure for data analyses:
     MyStat1 = stats_t( Fs1 );
     MyStat2 = stats_t( Fs2 );
     
    Hmu = figure(1);
    % set(Hmu,'WindowStyle','docked'); 
    % raw/(0.45*1e5)  [V]unit x [W/V] = [W]unit
   
%     [t1,Datas_mu1]   = MyStat1.average( raw/(0.45*1e5) );
%     Datas_std1  = MyStat1.standard_dev( raw/(0.45*1e5) ) ;   
%     [t2,Datas_mu2]   = MyStat2.average( raw'/(0.45*1e5) );
%     Datas_std2  = MyStat2.standard_dev( raw'/(0.45*1e5) );
       
    [Datas_mu1,Datas_std1,Datas_mu2,Datas_std2] = AverageDataBothWays( raw/(0.45*1e5) );
    t1 = (1/Fs2)*(1:length(Datas_mu1));
    t2 = (1/Fs1)*(1:length(Datas_mu2));
    
    
    % define plot colors
    
    h1 = subplot(2,2,1);
    set(h1,'parent',Hmu);
    plot_areaSTD(t1*1e3,1e6*Datas_mu1,1e6*sqrt(Ephoton*Fs1*Datas_mu1), h1,[243 169 114]./255,0.5);
    hold(h1,'on')
    plot_areaSTD(t1*1e3,1e6*Datas_mu1,1e6*Datas_std1, h1,[128 193 219]./255,0.5);
    ylabel('\mu (\mu W)over short ')
    xlabel('time (ms)')
    h2 = subplot(2,2,2);
    hold(h2,'on')
    plot_areaSTD(t2*1e6,1e6*Datas_mu2,1e6*sqrt(Ephoton*Fs1*Datas_mu2), h2,[243 169 114]./255,0.5);
    hold(h2,'on')
    plot_areaSTD(t2*1e6,1e6*Datas_mu2,1e6*Datas_std2, h2,[128 193 219]./255,0.5);
    ylabel('\mu (\mu W)over  long')
    xlabel('time (\mu s)')
    
% line(t1*1e3,1e6*Datas_mu1,'Color','r'); hold on 
%     ylabel('\mu (\mu W)over short ')
%     xlabel('time (ms)')
%     ax1 = gca; % current axes
%     set(ax1,'XColor','r');
%     set(ax1,'YColor','r');
%     ax1_pos = get(ax1,'Position'); % position of first axes
%     ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% line(t2*1e6,1e6*Datas_mu2,'Parent',ax2,'Color','k')
%     set(ax2,'Ylim',get(ax1,'Ylim'));
%     ylabel('\mu (\mu W)over  long')
%     xlabel('time (\mu s)')

    subplot(2,2,[3 4])
    % psdx unit: [raw^2/Hz^2] = [W^2/Hz^2] (cf Equation 3)
[freq1 , psdx1 ] = MyStat1.PowSpecDens( raw/(0.45*1e5) ) ;
[freq2 , psdx2 ] = MyStat2.PowSpecDens( raw'/(0.45*1e5) ) ;
% averaging over the many realizations:
psdx1 = mean(psdx1,2);
psdx2 = mean(psdx2,2);

%  figure;plot(freq1,10*log10(Fs1*psdx1)),
%  figure;plot(freq2,10*log10(Fs2*psdx2),'Color','r')
% e1  = MyStat1.Energy_t( raw/(0.45*1e5) );
% ee1 = MyStat1.Energy_psd( psdx1 ); 

line(freq2,10*log10(psdx2),'Color','r')
    ylim([-300 -80])
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    ax1 = gca; % current axes
    set(ax1,'XColor','r');
    set(ax1,'YColor','r');
ax1_pos = get(ax1,'Position'); % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(freq1*1e-6,10*log10(psdx1),'Parent',ax2,'Color','k')
     set(ax2,'Ylim',get(ax1,'Ylim'));
     set(ax2,'Xlim',[min(freq1*1e-6) max(freq1*1e-6)]);
    xlabel('Frequency (MHz)')
    ylabel('Power (dB)')

   
    