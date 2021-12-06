function D = ReadHighFinessFile(path)
%This function created by Maïmouna Bocoum 29-10-2021 reads Highfiness files

D = importdata(path);

Data = D.data ;

Header = D.textdata ; 

% fileID = fopen(path);
% 
% 
% 
% Data_out = fscanf(fileID,'%d');
% 
% 
% fclose(fileID);

end

