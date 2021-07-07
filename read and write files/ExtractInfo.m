function MyString_ = ExtractInfo(MyString,StringBefore,StringAfter)

% function intended to extract index in complicated file name. 
% using regexpi, here are the string interpretation:

% . Any single character

% run a loop to skip '*' seperated componant:
T_star = findstr(StringBefore,'*')
splitcells = regexp( StringBefore , '*','split')
% [L,idx,id0] = unique(splitcells)
% 
% Ip = findstr(MyString,StringBefore)
MyString_ = MyString ;

for n_loop = 1:length(splitcells)

    MyString_ = extractAfter( MyString_ ,splitcells{n_loop}) ;
      
end

% extract information until "StringAfter":

MyString_ = extractBefore( MyString_ , StringAfter ) ;

   


end

