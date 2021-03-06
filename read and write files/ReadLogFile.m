function [NB,nuX0,nuZ0] = ReadLogFile(SubFolderName)



% gets todays path:

% SubFolderName  = generateSubFolderName('D:\Data\Mai');
% SubFolderName = 'D:\Data\Mai\2020-01-14\NbXNbZ';
LogFile = importdata([SubFolderName,'\LogFile.csv']);

TypeOfSequence      = LogFile.colheaders{:,1};
Volt                = LogFile.colheaders{:,2} ;
FreqSonde           = LogFile.colheaders{:,3};
NbHemicycle         = LogFile.colheaders{:,4};
Tau_cam             = LogFile.colheaders{:,5};
DurationWaveform    = LogFile.colheaders{:,6} ;
Prof                = LogFile.colheaders{:,7};
Nevent              = LogFile.colheaders{:,8} ;
NTrig               = LogFile.colheaders{:,9} ;
X0                  = LogFile.colheaders{:,10};
nuX0                = str2double( LogFile.colheaders{:,11} ); 
nuZ0                = str2double( LogFile.colheaders{:,12} ) ;

NB = LogFile.data(:,13:16) ;




end

