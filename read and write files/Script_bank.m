%% script Bank

%% read multiple files
% [Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
% Nfiles_ = length(Filename);
% if Nfiles_==1
%    Filename = {Filename};
% end
% 
% %% test ReadDatFile as ascii
% I = ReadDatFile( [Foldername,Filename] )

%% read banary file %%
figure(3)
%Data_zero = ReadBinFile( 'D:\Data\Mai\2020-09-02\Caxis_CW1169mA_noPump_525mW_zero_12_02_07.dat' ,'matrix_DBL');
Data = ReadBinFile( 'D:\Data\Mai\2020-08-27\PD_state3_500mW_pump100W_5_04_53.dat' ,'matrix_DBL');
%Cal = 0.290/0.0117 ; % W/Volt-2020-07-17
%Cal = 0.347/0.0523 ; % W/Volt-2020-08-27   - PD state 1
Cal = 0.4/0.0141 ;   % W/Volt-2020-08-27   - PD state 2
% density 09-02-2020 : 0.0584
% Cal without density:
%Cal = 0.525/0.0916 ; % W/Volt-2020-09-02
hold on
time = 1e6*(1:size(Data,2))/(1e7); % time axis in us
plot( time - mean(time) + 59.15 ,Cal*Data(1,:) - mean(Cal*Data(1,1:30)) )
% xlim([-100 300])
title('Pump 100 W - influence delais injection')
xlabel('Time (\mu s)')
ylabel('Output Power (W)')

% legend('pump(W)','seed(W)','simu amplified(W)','exp amplified(W)')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
%%

%  seed = [0.280,0.5,0.6,0.7,0.8,0.9,1,1.3];
%  ouput_manip = [3.995,5.1,5.526,5.871,6.5,6.691,7.32,8.203];
% % ouput_simu = [2.918,3.23,3.361,3.4,3.617,3.722,3.837,4.169];
% ouput_simu = [5.292,5.88,6.127,6.363,6.59,6.803,7.025,7.659];
% figure; plot(seed,ouput_manip) ; hold on ;  plot(seed,ouput_simu) ;
% legend('manip(W)','simu(W)')
% xlabel('seed power')
% ylabel('ouput power (plateau)')

%% differential 

%%






