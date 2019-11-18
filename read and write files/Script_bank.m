%% script Bank

%% read multiple files
[Filename,Foldername] = uigetfile('*.tiff','MultiSelect','on');
Nfiles_ = length(Filename);
if Nfiles_==1
   Filename = {Filename};
end