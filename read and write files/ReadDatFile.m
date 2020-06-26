function [ Data_out ] = ReadDatFile( path )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fileID = fopen(path);

% change into line by line regular expression

Data_out = fscanf(fileID,'%d');

fclose(fileID);

end

