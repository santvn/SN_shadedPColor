function [ax,h_c,h_s,cx]=SN_shadedPColor(varargin)
%SN_SHADEDPCOLOR creates shaded pcolor plot using one dataset as a color and
%   one dataset as a shading on top
%
%   SN_SHADEDPCOLOR(C,S) uses subplots to create the main pcolor plots and
%   colorbars. C is a matrix of data for color plots. S is a matrix used
%   for shading
% 
%   SN_SHADEDPCOLOR(X,Y,C,S) uses X and Y (vectors or matrices) to
%   indicates the x- and y-coordinates for each point in matrices C and S.
%   
%   [AX,HC,HS,CX] = SN_SHADEDPCOLOR(...) returns the graphic handles for
%   the major graphical axes AX, color plot HC, shaded plot HS, and
%   colorbars for color and shading as CX.
%   
%   SN_SHADEDPCOLOR(C,S,'OPTION',OPTION_VAL,...) allows the users to
%   specify options for customization of the shaded plot
%   The OPTIONS are:
%       'cMin' :  minimum value for color plot
%       'cMax' :  maximum value for color plot
%       'cLim' or 'cRange' or 'cAxis' :  minimum & maximum values for color plot
%       'sMin' :  minimum value for shaded plot
%       'sMax' :  maximum value for shaded plot
%       'sLim' or 'sRange' or 'sAxis' :  minimum & maximum values for shaded plot
%       'cmap' or 'colormap':  specific color map for the color plot
%       'shadingExp': shading exponent for shading, this provides a sharp or
%           smooth curve to the shading
%       'shadingFact': shading factor for shading to make the shading light
%           or dark
%       'reverseShade': indicator for inverse shading
%
% See also PCOLOR
%
% Created by San Nguyen 2015 09 16 (please report any known bugs to
%                                   stn004@ucsd.edu)
%

if nargin<2
    error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
end
    
persistent argsNameToCheck;
if isempty(argsNameToCheck);
    argsNameToCheck = {...
        'cMin',... %1
        'cMax',... %2
        'cRange',... %3
        'cLim',... %4
        'cAxis',... %5
        'sMin',... %6
        'sMax',... %7
        'sRange',... %8
        'sLim',... %9
        'sAxis',... %10
        'cmap',... %11
        'colormap',... %12
        'shadingExp',... %13
        'shadingFact',... %14
        'reverseShade'... %15
        };
                        
end

shading_exp = 3;
shading_fac = 0.5;
cmap = get(0,'defaultFigureColormap');
c_min = NaN;
c_max = NaN;
s_min = NaN;
s_max = NaN;
reverseShade = false;

n = 4;
if nargin<4
    n = nargin;
end

n_valid_data = 0;
for i = 1:n
    if ischar(varargin{i}) && isvector(varargin{i})
        break;
    end
    n_valid_data = n_valid_data + 1;
    
end
if ismember(n_valid_data,[0 1 3])
    error('MATLAB:SN_shadedPColor:missingArgs',...
        'Missing input arguments');
end

switch n_valid_data
    case 2
        c = varargin{1};
        s = varargin{2};
        x = 1:size(c,2);
        y = 1:size(c,1);
    case 4
        x = varargin{1};
        y = varargin{2};
        c = varargin{3};
        s = varargin{4};
end

if size(c,1) ~= size(s,1) || size(c,2) ~= size(s,2) ||...
        numel(x)~=size(c,2) || numel(y) ~= size(c,1)
    error('MATLAB:SN_shadedPColor:invalidInputs',...
        'Data must have the same dimensions');
end

index = n_valid_data+1;
n_items = nargin-n_valid_data;
while (n_items > 0)
    argsMatch = strcmpi(varargin{index},argsNameToCheck);
    i = find(argsMatch,1);
    if isempty(i)
        index = index +1;
        n_items = n_items-1;
        continue;
    end
    
    switch i
        case 1 % cMin
            if n_items == 1
                error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
            end
            c_min = varargin{index+1};
            if isnan(c_min) || isinf(c_min) || ~isscalar(c_min)
                error('MATLAB:SN_shadedPColor:invCMin',...
                    'Invalid value for cMin');
            end
            index = index +2;
            n_items = n_items-2;
        case 2 % cMax
            if n_items == 1
                error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
            end
            c_max = varargin{index+1};
            if isnan(c_max) || isinf(c_max) || ~isscalar(c_max)
                error('MATLAB:SN_shadedPColor:invCMax',...
                    'Invalid value for cMax');
            end
            index = index +2;
            n_items = n_items-2;
        case {3,4,5} % cRange
            if n_items == 1
                error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
            end
            
            if sum(isnan(varargin{index+1})) ||...
                    sum(isinf(varargin{index+1})) ||...
                    numel(varargin{index+1})~=2
                error('MATLAB:SN_shadedPColor:invCRange',...
                    'Invalid value for cRange');
            end
            
            c_min = varargin{index+1}(1);
            c_max = varargin{index+1}(2);
            
            index = index +2;
            n_items = n_items-2;
        case 6 % sMin
            if n_items == 1
                error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
            end
            s_min = varargin{index+1};
            if isnan(s_min) || isinf(s_min) || ~isscalar(s_min)
                error('MATLAB:SN_shadedPColor:invSMin',...
                    'Invalid value for sMin');
            end
            index = index +2;
            n_items = n_items-2;
        case 7 % sMax
            if n_items == 1
                error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
            end
            s_max = varargin{index+1};
            if isnan(s_max) || isinf(s_max)|| ~isscalar(s_max)
                error('MATLAB:SN_shadedPColor:invSMax',...
                    'Invalid value for SMax');
            end
            index = index +2;
            n_items = n_items-2;
        case {8,9,10} % sRange
            if n_items == 1
                error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
            end
            
            if sum(isnan(varargin{index+1})) ||...
                    sum(isinf(varargin{index+1})) ||...
                    numel(varargin{index+1})~=2
                error('MATLAB:SN_shadedPColor:invSRange',...
                    'Invalid value for sRange');
            end
            
            s_min = varargin{index+1}(1);
            s_max = varargin{index+1}(2);
            
            index = index +2;
            n_items = n_items-2;
        case {11,12} % cmap
            if n_items == 1
                error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
            end
            cmap = varargin{index+1};
            if sum(isnan(cmap(:))) || size(cmap,2)~=3 || sum(cmap(:)>1) || sum(cmap(:)<0)
                error('MATLAB:SN_shadedPColor:cMap',...
                    'Invalid colormap');
            end
            index = index +2;
            n_items = n_items-2;
        case 13 % shading_exp
            if n_items == 1
                error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
            end
            shading_exp = varargin{index+1};
            if ~(isscalar(shading_exp) && isnumeric(shading_exp)) ||...
                    isnan(shading_exp) || isinf(shading_exp)
                error('MATLAB:SN_shadedPColor:shadingExp',...
                    'Shading exponent should be finite and greater than 0');
            end
            index = index +2;
            n_items = n_items-2;
        case 14 % shadingFact
            if n_items == 1
                error('MATLAB:SN_shadedPColor:missingArgs',...
                    'Missing input arguments');
            end
            shading_fac = varargin{index+1};
            if ~(isscalar(shading_fac) && isnumeric(shading_fac)) ||...
                    isnan(shading_fac) || isinf(shading_fac)
                error('MATLAB:SN_shadedPColor:shadingFact',...
                    'Shading factor should be finite and greater than 0');
            end
            index = index +2;
            n_items = n_items-2;
        case 15 % reverseShade
            reverseShade = true;
            index = index +1;
            n_items = n_items-1;
    end
end

if ~isnan(c_min) && ~isnan(c_max) && c_min==c_max
    error('MATLAB:SN_shadedPColor:invCLim',...
        'cMin has to be less than cMax');
end

if isnan(c_min) 
    c_min = nanmin(c(:));
end

if isnan(c_min) 
    c_min = 0;
end

if isnan(c_max)
    c_max = nanmax(c(:));
end
if isnan(c_max)
    c_max = 1;
end

if c_min > c_max
    tmp = c_min;
    c_min = c_max;
    c_max = tmp;
end
% keyboard;
if ~isnan(s_min) && ~isnan(s_max) && s_min==s_max
    error('MATLAB:SN_shadedPColor:invSLim',...
        'sMin has to be less than sMax');
end

if isnan(s_min) 
    s_min = nanmin(s(:));
end

if isnan(s_min) 
    s_min = 0;
end

if isnan(s_max)
    s_max = nanmax(s(:));
end
if isnan(s_max)
    s_max = 1;
end

if s_min > s_max
    tmp = s_min;
    s_min = s_max;
    s_max = tmp;
end

c = (c-c_min)/(c_max-c_min);
c(c<0) = 0;
c(c>1) = 1;

s = (s-s_min)/(s_max-s_min);
s(s<1/size(get(0,'defaultfigurecolormap'),1)) = 1/size(get(0,'defaultfigurecolormap'),1);
s(s>1) = 1;
s = -s;

ax = subplot('position',[0.075 0.15 .785 .8]);
h_c = pcolor(x,y, c);

caxis([-1 1]);
if reverseShade
    colormap([flipud(gray(size(cmap,1)));cmap]);
else
    colormap([(gray(size(cmap,1)));cmap]);
end


hold on;
h_s = pcolor(x,y, s);

alpha(h_s,shading_fac*(-s).^(shading_exp));
set(h_c,'FaceColor','flat');
set(h_s,'FaceColor','flat');
set(h_s,'FaceAlpha','flat');
% ylabel('depth [m]');
% xlabel('time [yearday 2015, UTC]');


cx(1) = subplot('position',[0.87 0.15 .01 .8]);
x = [0; 1];
y = linspace(c_min,c_max,floor(size(colormap,1)/2));
c = repmat(linspace(0,1,floor(size(colormap,1)/2)),[2 1])';
contourf(x,y,c,...
    floor(size(colormap,1)/2),'linestyle','none');
caxis([-1 1])
set(gca,'yAxisLocation','right','tickLength',[0.01 0.005],'tickDir','in','xtick',[])
% ylabel('temperature [$^\circ$C]');

cx(2) = subplot('position',[0.94 0.15 .01 .8]);
x = [0; 1];
y = linspace(s_min,s_max,floor(size(colormap,1)/2));
c = repmat(linspace(0,1,floor(size(colormap,1)/2)),[2 1])';
c = -(c.^(shading_exp+1));
contourf(x,y,c,...
    floor(size(colormap,1)/2),'linestyle','none');
caxis([-1 1])
set(gca,'yAxisLocation','right','tickLength',[0.01 0.005],'tickDir','in','xtick',[])
% ylabel('$\partial \rho/\partial z$ [kg/m$^4$] $\times 10^{-3}$');
axis ij;
grid on;
box on;
end