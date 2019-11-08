%%% main %%%
clearvars

%% read Log - file
pathLog = 'D:\Data\Mai\2018-11-22\Scans\LogScanSansBascule_1.dat';
LOG = importdata(pathLog);
column = 1 ;
param = LOG.data(:,column);

%% read Data
pathData = 'D:\Data\Mai\2018-11-22\Scans\ScanSansBasucle_1.dat';
[Header,Data] = ReadDataFile(pathData);
D = importdata(pathData);

% gage header
Fs  =   str2double( regexp(Header{3},['\d+'],'match') );
N   =   str2double( regexp(Header{7},['\d+'],'match') );
Nav =   str2double( regexp(Header{11},['\d+'],'match') );

t = (1:N)*(1/Fs);

figure; imagesc(t*1e6,param,Data)
ylabel(LOG.colheaders{column})
xlabel('time (\mu s)')


