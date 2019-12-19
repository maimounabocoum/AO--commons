% ----------------------------------------------------------------------- %
% Function plot_areaerrorbar plots the mean and standard deviation of a   %
% set of data filling the space between the positive and negative mean    %
% error using a semi-transparent background, completely customizable.     %
%                                                                         %
%   Input parameters:                                                     %
%       - data:     Data matrix, with rows corresponding to observations  %
%                   and columns to samples.                               %
%       - options:  (Optional) Struct that contains the customized params.%
%           * options.handle:       Figure handle to plot the result.     %
%           * options.color_area:   RGB color of the filled area.         %
%           * options.color_line:   RGB color of the mean line.           %
%           * options.alpha:        Alpha value for transparency.         %
%           * options.line_width:   Mean line width.                      %
%           * options.x_axis:       X time vector.                        %
%           * options.error:        Type of error to plot (+/-).          %
%                   if 'std',       one standard deviation;               %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval.              %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       data = repmat(sin(1:0.01:2*pi),100,1);                            %
%       data = data + randn(size(data));                                  %
%       plot_areaerrorbar(data);                                          %
% ----------------------------------------------------------------------- %
%   Author:  Maïmouna Bocoum                                      %
%   Date:    12/12/2019                                                   %
%   E-mail:  maimmouna.bocoum (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
function [] = plot_areaSTD(x,data,sigma,Ha,c,alpha)

% Ha : target axes to plot

    % Default options
 

        options.color_area = c ;    % Blue theme
        options.color_line = [ 52 148 186]./255;
        %options.color_area = [243 169 114]./255;    % Orange theme
        %options.color_line = [236 112  22]./255;
        options.alpha      = alpha;
        options.line_width = 2;
        options.error      = 'std';


    options.x_axis = x(:);
    
    % Computing the mean and standard deviation of the data matrix
 
    % Plotting the result
    
    x_vector = [options.x_axis', fliplr(options.x_axis')];
    patch = fill(x_vector, [data+sigma,fliplr(data-sigma)], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    set(patch,'parent',Ha)
    
    hold on;
    plot(options.x_axis, data, 'color', options.color_line, ...
        'LineWidth', options.line_width,'parent',Ha);
    hold off;
    
end
