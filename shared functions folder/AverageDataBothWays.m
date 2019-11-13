function [Datas_mu1,Datas_std1, Datas_mu2 , Datas_std2] = AverageDataBothWays( raw )
% To maintain the default normalization 
% while specifying the dimension of operation, set w = 0 in the second argument
% of var function

Datas_mu1 = mean(raw,1);
Datas_std1 = sqrt(var(raw,0,1));


Datas_mu2 = mean(raw,2);
Datas_std2 = sqrt(var(raw,0,2));


end

