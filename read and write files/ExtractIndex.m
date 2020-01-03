function Index = ExtractIndex(MyString)



%% typical nume blabkla_1234.tiff
% extractinf : 1234 as double

Is = findstr(MyString,'_');

Ip = findstr(MyString,'.tif');

Index = str2double( MyString((Is(end)+1):(Ip(end)-1)) );









end

