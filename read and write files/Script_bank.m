%% script Bank

%% read multiple files
[Filename,Foldername] = uigetfile('D:\Data\Mai\2020-09-10\*.dat','MultiSelect','on');
Nfiles_ = length(Filename);
if Nfiles_==1
   Filename = {Filename};
end


% test ReadDatFile as ascii
% Data_zero   = ReadBinFile( [Foldername,Filename{2}],'matrix_DBL');
% Data        = ReadBinFile( [Foldername,Filename{1}],'matrix_DBL');

%% read banary file %%
cc = jet(Nfiles_);
for loop = 1:Nfiles_
 figure(1)
 %Data_zero   = ReadBinFile( 'D:\Data\Mai\2020-11-19\Calibration_zeroCW_5_42_43.dat' ,'matrix_DBL');
% Data        = ReadBinFile( 'D:\Data\Mai\2020-11-23\time response using AM main modulation\AMmain10kHz_pump100W_4_54_37.dat' ,'matrix_DBL');

Data        = ReadBinFile( [Foldername,Filename{loop}],'matrix_DBL');
Results(loop,:) = Data(1,:);
% figure(2)
% plot((time(1500:end)-500)/4,medfilt1(Results2',5,[],1))
% xlabel('time (\mu s)')
% ylabel('Power W)')
% cb =  colorbar ;
% ylabel(cb,'normalized relative power')
% title('Influence of frequency mismatch')
% set(findall(gcf,'-property','FontSize'),'FontSize',28)



%Cal = 0.290/0.0117 ; % W/Volt-2020-07-17
%Cal = 0.347/0.0523 ; % W/Volt-2020-08-27   - PD state 1
%Cal = 0.4/0.0141 ;   % W/Volt-2020-08-27   - PD state 2
Cal=1; % AsGa type II - 2020-09-07 
% density 09-02-2020 : 0.0584
% Cal without density:
%Cal = 0.525/0.0916 ; % W/Volt-2020-09-02
Cal = 1;
hold on
time = 1e6*(1:size(Data,2))/(20000); % time axis in us
[a,b] =max(Data(1,1400:1800))
plot( time(200:2000) , Data(1,200:2000)-Data(1,1000) , 'color',cc(loop,:),'linewidth',1)
%plot( Data(1,1500:1800), 'color',cc(loop,:),'linewidth',3)
title('Negative Gain')
xlabel('Time (\mu s)')
ylabel('Volt')

% legend('pump(W)','seed(W)','simu amplified(W)','exp amplified(W)')
end
set(findall(gcf,'-property','FontSize'),'FontSize',14)
%%
figure(2);
hold off
y2 = mean( Results(8000:14000,:) , 1 ) ; 
y1 = mean( Results(2000:4000,:) , 1 ) ;
df = [0:10,20:10:100,200,300,400];
figure(2); plot( 1e3./df , y1 ,'linewidth', 6)
xlabel('\tau (\mu s)')
ylabel('signal photodiode (a.u)')
c1 = min(find(y1>=( max(y1) + min(y1) )/2 ));
hold on ; plot( 1e3./df , y2 ,'linewidth', 6)
xlabel(' (\Delta f)^{-1} (\mu s)')
ylabel('signal photodiode (a.u)')
c2 = min(find( y2>=( max(y2) + min(y2) )/2 ) );
grid on
legend('pump 0W - \tau_{1/2} = 200 \mu s','pump 100W - \tau_{1/2} = 12.5 \mu s')
title('Influence of pump power on crystal response time')
set(findall(gcf,'-property','FontSize'),'FontSize',22)

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






