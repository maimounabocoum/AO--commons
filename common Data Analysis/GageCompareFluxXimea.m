%% gage loop analysis
%clearvars
%load('Q:\datas\2019-12-10\Ref_100Hz_41uW_filter_main2uW_RepHz_100_16h27_53.mat');
% load('D:\Data\Mai\2020-03-06\6MHz-8 Volt\Water_Fus_6_Volt_8_18h22_47');
%% =============== add subfolders ===============  %%

%  addpath('..\..\AO--commons\Plottings')
%  addpath('..\..\AO--commons\shared functions folder')
%  addpath('gui');
%  addpath('sequences');
%  addpath('..\..\AO--commonssubfunctions');

%% please run Ximea analysis prior to this script

%% =============== screen out datas ===============  %%

     load([Foldername,Filename{nfile}])
     h = 6.6e-34;
     lambda = 780e-9;
     Ephoton = h*(3e8/lambda);
     IntDoor = zeros(1,length(t2));
     ts = 235e-6;
     tr = ts + 50e-6;
     Iextract_s = find(t2 >= ts & t2 <= ts + MyXimea.IntegrationTime);
     Iextract_r = find(t2 >= tr & t2 <= tr + MyXimea.IntegrationTime);
     
     % extraction of both signal and reference
     SignalCamera = raw(Iextract_s,:);
     SignalRef   = raw(Iextract_r,:);
     
     % just to screening we define the door on signal and reference
     IntDoor(Iextract_s) = 1;
     IntDoor(Iextract_r) = 1;
     
    % check data extraction using door 
    Hmu = figure(1);
    set(Hmu,'WindowStyle','docked'); 
    plot(t2*1e6,mean( raw ,2)); hold on ;
    plot(t2*1e6,IntDoor); hold off ;
    xlim([200 300])
    % plot(t2(Iextract)*1e6,mean(SignalCamera ,2))
  %%  evaluation of signal ratio
  
 s_pd = ( sum(SignalRef,1) - sum(SignalCamera,1) )./sum(SignalRef,1) ;
 

    
%% plot Ptot
    figure(3) ; hold on
    plot(repmat(P(nfile),1,length(s_pd(4:end))),100*s_pd(4:end),'o')
    xlabel('shot index')
    title('measured power')
    % legend('Holograhie','Photodiode') ; 
    
    
    %%
   

 %% Bessel functions
 phi = 0:0.01:2*pi;
    Hmu = figure(2);
    set(Hmu,'WindowStyle','docked');
    bsum = 100*besselj(0 ,phi).^2;
    for loop = 0:4
    plot( phi , 100*besselj(loop ,phi).^2 ); hold on
    if loop>0
    bsum = bsum + 2*100*besselj(loop ,phi).^2;
    end
    end
    plot(phi,bsum); hold off;
    
    legend('J_0','J_1','J_2','J_3','J_4')
    xlabel('\phi(Rad)')
    ylabel('efficiency [%]')
    title('J_{n}(\phi)^2')
    hold off








    