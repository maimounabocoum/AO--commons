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
Data = ReadBinFile( 'D:\Data\Mai\2020-06-24\test.dat' ,'matrix_DBL')
%%