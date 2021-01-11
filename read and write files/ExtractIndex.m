function Index = ExtractIndex(MyString,mode)



%% typical nume blabkla_1234.tiff
% extractinf : 1234 as double

switch mode
    case 'labview'
Is = findstr(MyString,'_');

Ip = findstr(MyString,'.tif');

Index = str2double( MyString((Is(end)+1):(Ip(end)-1)) );
    case 'ximea'
Ip = findstr(MyString,'.tif');

Index = str2double( MyString(1:(Ip(end)-1)) );        
end









end

