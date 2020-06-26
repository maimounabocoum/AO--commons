function [ Data ] = ReadBinFile( File , type )
%created by maimouna bocoum 6-25-2020

fileID = fopen(File);
Data = [];

switch type
    
    case 'matrix_DBL'
        while ~feof(fileID)
        % retreive matrix dimension
        s = fread(fileID,2,'int','ieee-be') ; 
        % skip matrix dimension of type double
        Data = fread(fileID,'double','ieee-be');
        
        % resize data
        Data = reshape(Data',s(2),s(1))' ;
        end
        
    case 'vector_DBL'
        while ~feof(fileID)
        s = fread(fileID,1,'int','ieee-be') ;
        Data = fread(fileID,'double','ieee-be');
        end
end




fclose(fileID);

end

