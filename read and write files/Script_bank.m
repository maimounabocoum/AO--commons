%% script Bank

%% read multiple files
[Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
Nfiles_ = length(Filename);
if Nfiles_==1
   Filename = {Filename};
end

%% test ReadDatFile as ascii
I = ReadDatFile( [Foldername,Filename] )

%% read banary file %%
Data = ReadBinFile( 'D:\Data\Mai\2020-06-26\mydata_12_02_45.dat' ,'matrix_DBL')
%%