function varargout = PlotDistBands(y,varargin)

% PlotDistBands
%
% Plots median and percentile bands for a matrix.
%
% Usage: 
%   PlotDistBands(y)
%   PlotDistBands(x,y)
%   PlotDistBands(...,OptionsStructure,...)
%   PlotDistBands(...,'PropertyName',PropertyValue,...)
%   h = PlotDistBands(...)
%
%
% [required input]:
%   - y
%   Data matrix. Percentiles will be computed along its first dimension.
%                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                     
% [Options]:
%
%   - Bands2Show
%   Percent intervals to be shown, centered around median.
%   Default: [50, 60, 70, 80, 90]
%
%   - MedianColor
%   Color of median line.
%   Default: [0, 0, 0.7]
%
%   - ShadeColor
%   Base color for bands.
%   Default: [0.2, 0.6, 0.5]
%
%   - LineWidth
%   Width of median line.
%   Default: 1.5
%
%   - isZeroLine
%   If 1, plots the zero line. If 0, it does not plot a zero line.
%   Default: 1
%
%   - ZeroLineColor
%   Color for the zero line.
%   Default: 'k'
%
%   - ZeroLineStyle
%   Style for the zero line.
%   Default: '-'
%
%   - ZeroLineWidth
%   Width for the zero line.
%   Default: 0.5
%
%   - tid
%   x axis values.
%   Default: 1:T
%
% See also:
% AdvPlot, ColorPallette
%


%% ------------------------------------------------------------------------

%% 1. Check main inputs
if nargin == 0
    error('Error: No arguments are specified.');
elseif nargin == 1 % PlotDistBands(y)
    o = struct;
else
    if isnumeric(varargin{1}) % PlotDistBands(x,y, ...)
        x = y;
        y = varargin{1};
        varargin(1) = [];
    end
    if isstruct(varargin{1}) % PlotDistBands(y, o, ...)
        o = varargin{1};
        varargin(1) = [];
    end
    for jo = 1:(length(varargin)/2) % Load Options
        o.(varargin{(jo-1)*2+1}) = varargin{jo*2};
    end
end

if ~exist('x','var')
    x = 1:size(y,2);
end

%% 2. Apply default options to unspecified ones
oDefault.Bands2Show           = [50,60,70,80,90];
oDefault.AltData              = [];
oDefault.ShowLegend           = 0;
oDefault.LegendLocation       = 'Best';
oDefault.LegendOrientation    = 'vertical'; %'vertical','horizontal'
oDefault.LegendWithBands      = 0;
oDefault.LineColor            = ColorPallette;
oDefault.ShadeColor           = [0.72, 0.77, 0.82];
oDefault.ShadeColorBrightness = 0.9;
oDefault.ShadeFactors         = [0.1, 0.65]; % shade factors at 50 and 90%
oDefault.isZeroLine           = 1;
oDefault.ZeroLineColor        = 'k';
oDefault.ZeroLineStyle        = '-';
oDefault.ZeroLineWidth        = 0.5;

% MedianColor = [0,0,0.7];
% ShadeColor = [0.2,0.6,0.5];

oList = fieldnames(oDefault);
for jo = 1:length(oList)
    oName = oList{jo};
    if ~isfield(o,oName)
        o.(oName) = oDefault.(oName);
    end
end


%% 3. Additional options, if not set before
% o.LineColor = [0.3,0.75,0.3;o.LineColor];
if ~isfield(o,'ShadeColor')
    if ~isfield(o,'ShadeColorBase')
  %   o.ShadeColorBase = [0.70,0.70,0.7]; % grey
  %   o.ShadeColorBase = [0.70,0.80,0.85]; % light blue 
  %   o.ShadeColorBase = [0.60,0.70,0.75]; % grey/Blue
        o.ShadeColorBase = [0.15,0.25,0.75]; % Blue
    end
    if ~isfield(o,'ShadeCLWeight')
        o.ShadeCLWeight = 0.0;
    end
%   o.ShadeColor = [0.65,0.75,0.8]*0.95+0.05*o.LineColor(1,:);
    o.ShadeColor = o.ShadeColorBrightness * o.ShadeColorBase * ( 1 - ...
                   o.ShadeCLWeight) + o.ShadeCLWeight*o.LineColor(1,:);
end

%% 4. Plot zero line
if o.isZeroLine
    h.ZeroLine = plot(x, zeros(1,size(y,2)), o.ZeroLineStyle,...
                      'Color'     ,o.ZeroLineColor,...
                      'LineWidth' ,o.ZeroLineWidth);
    hold on
end

%% 5. Plot bands
o.Bands2Show   = sort(o.Bands2Show,'descend');
nBands         = length(o.Bands2Show);
InitHold       = ishold;
BandsData      = zeros(nBands*2, size(y,2));
BandColorSlope = [-1,1] * o.ShadeFactors' / ([1,-1] * o.Bands2Show([1,nBands])');
BandColorCt    = o.ShadeFactors(1) - BandColorSlope * o.Bands2Show(nBands);
for jB = 1:nBands
    Band     = o.Bands2Show(jB);
    BandPath = prctile(y, 50 + Band/2 * [-1,+1]); % Upper band and lower band
    BandColor = o.ShadeColor +...
               (1 - o.ShadeColor)*(BandColorCt + BandColorSlope * Band);
    h.Bands(jB) = fill([x,x(end:-1:1)], [BandPath(1,:),BandPath(2,end:-1:1)],...
                       BandColor, 'EdgeColor', BandColor); 
    hold on
    BandsData((jB-1)*2+[1,2], :) = BandPath;
end
h.BandsData = BandsData;


%% 6. Plot Median and alternative data
%  AdvPlot.m is required

YData = [prctile(y,50); o.AltData];
if o.ShowLegend
    if ~isfield(o,'LegendString')
        o.LegendString{1} = 'Median';
        for j = 2:size(YData,1)
            o.LegendString{j} = sprintf('Alt %.0f',j-1);
        end
        for j = 1:nBands*o.LegendWithBands
            o.LegendString{size(YData,1)+j} = sprintf('%.0f%%',o.Bands2Show(j));
        end
    end
    if o.LegendWithBands
        o.LegendItems = h.Bands;
    end
end

h.Lines = AdvPlot(x, YData, o);

o.ShowLight = h.Lines.Options.ShowLight;
h.LegendObj = h.Lines.LegendObj;
if o.ShowLegend && o.ShowLight
    h.LegendLines = h.Lines.LegendLines;
end


%% 7. Postamble
if ~InitHold
    hold off
end
xlim([x(1),x(end)])

h.XData   = x;
h.YData   = YData;
h.Options = o;
if nargout == 1,varargout = {h};end

%% ------------------------------------------------------------------------
