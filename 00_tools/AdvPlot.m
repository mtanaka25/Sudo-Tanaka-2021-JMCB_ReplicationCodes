function varargout = AdvPlot(y,varargin)
% AdvPlot
%
% Plots several lines with adavanced opitions, such as showing a zero line
% and gradating graphes.
% 
% Syntax:
%   AdvPlot(y)
%   AdvPlot(x,y)
%   AdvPlot(...,OptionsStructure,...)
%   AdvPlot(...,'PropertyName',PropertyValue,...)
%   h = AdvPlot(...)

%% ------------------------------------------------------------------------

%% 1. Check main input
if nargin == 0
    error('Error: No arguments are specified.');
elseif nargin == 1  % AdvPlot(y)
    x = 1:size(y,2);
    o = struct;
else
    if isnumeric(varargin{1}) % AdvPlot(x,y, ...)
        x = y;
        y = varargin{1};
        varargin(1) = [];
    end
    if isstruct(varargin{1}) % AdvPlot((x,)y, o, ...)
        o = varargin{1};
        varargin(1) = [];
    end
    for jo = 1:(length(varargin)/2)  % AdvPlot((x,)y,(o,) 'Option', value, ...)
        o.(varargin{(jo-1)*2+1}) = varargin{jo*2};
    end
end

if ~exist('x','var')
  x = 1:size(y,2);
end


%% 2. Apply default options to unspecified ones
%oDefault.LineColor      = ColorPallette;
oDefault.LineColor      = [[ 255,   0,   0];... % Red
                           [   0,   0, 255];... % Blue
                           [ 255,   0,   0];... % Red
                           [   0,   0, 255];... % Blue
                           [   0,   0,   0];... % Black
                           [   0, 255,   0] ... % Green
                           ]/255; % normalize to unity
oDefault.LineStyle      = {'-', '-', '--', '--', ':', ':'};
oDefault.LineMarker     = {'none','none','none','none','none','none'};
oDefault.LineMarkerSize = [1, 1, 1, 3, 3, 3];
oDefault.LineWidth      = [2, 1, 1, 1, 1, 1];
oDefault.ShowLight      = 1;
oDefault.LightLayers    = 2;
oDefault.Brightness     = 0.3;
oDefault.ShowZeroLine   = 1;
oDefault.ZeroLineColor  = 'k';
oDefault.ZeroLineStyle  = '-';
oDefault.ZeroLineWidth  = 0.5;
oDefault.ShowLegend     = (size(y,1)>1);
oDefault.LegendLocation = 'Best'; 
oDefault.LegendOrientation = 'vertical'; %'vertical','horizontal'

oDefault.MarkerTicks    = max(1,round(40/size(y,2)));
% --- If there are not enough ticks for the lines with marker only then
%     the code automatically interpolates the data to add more marker ticks

oList = fieldnames(oDefault);
for jo = 1:length(oList)
    oName = oList{jo};
    if ~isfield(o,oName)
        o.(oName) = oDefault.(oName);
    end
end


%% 3. Set some auxiliary values
InitHold = ishold;
h.Lines  = zeros(1 + o.LightLayers * o.ShowLight, size(y,1));
xx       = x(1):1/o.MarkerTicks:x(size(y,2));

%% 4. Plot zero line
if o.ShowZeroLine
  h.ZeroLine = plot(xx,zeros(1,length(xx)),o.ZeroLineStyle,...
                    'Color'    , o.ZeroLineColor,...
                    'LineWidth', o.ZeroLineWidth);
  set(get(get(h.ZeroLine,'Annotation'),'LegendInformation'),...
      'IconDisplayStyle','off'); % Do not include ZeroLone in the legend
  if ~ishold
      hold on 
  end
else
  h.ZeroLine = [];
end


%% 5. Plot Line(s) with light effect
nLayers = 1 + o.ShowLight * o.LightLayers;

for j = 1 : size(y,1)
    yj      = y(j,:);
    yDiff   = diff(yj);
    yy      = zeros(1,length(xx));
    idx     = find(ismember(xx,x));
    yy(idx) = yj;
    
    for jM = 1:(o.MarkerTicks - 1)
        idx     = idx(1:size(y,2)-1) + 1;
        yy(idx) = yj(1:size(y,2)-1) + yDiff * jM / o.MarkerTicks;
    end
    for jLayer = 1 : nLayers
        if nLayers > 1
            jColor = o.LineColor(j,:) + ( 1 - o.LineColor(j,:) ) *...
                     o.Brightness * (jLayer - 1)/(nLayers - 1);
        else
            jColor = o.LineColor(j,:);
        end
        jSize = 1 - (jLayer - 1)/nLayers;
        h.Lines(jLayer,j) = plot(xx, yy,...
                                 'LineStyle'      ,o.LineStyle{j},...
                                 'Color'          ,jColor,...
                                 'MarkerFaceColor',jColor,...
                                 'Marker'         ,o.LineMarker{j},...
                                 'MarkerSize'     ,o.LineMarkerSize(j)*jSize,...
                                 'LineWidth'      ,o.LineWidth(j)*jSize);
        if ~ishold
            hold on
        end
    end
end

if isfield(o,'FontSize')
  set(gca,'FontSize',o.FontSize)
end

%% 6. Show Legend
if o.ShowLegend
    if ~isfield(o,'LegendString')
        for j = 1:size(y,1)
            o.LegendString{j} = sprintf('Line %.0f',j);
        end
    end
    h.LegendItems = h.Lines(1,:); % size(h,LegendItems) = size(y,1)
    if isfield(o,'LegendItems')
        h.LegendItems = [h.LegendItems,o.LegendItems];
    end
    nLegendItems = length(h.LegendItems);
    [h.Legend,h.LegendObj,~,~] = legend(h.LegendItems, o.LegendString,...
                                        'Location'   , o.LegendLocation,...
                                        'Orientation', o.LegendOrientation);
    for j = 1 : size(y,1) * o.ShowLight
        idx = nLegendItems + (j - 1) * 2 + 1;
        xx  = get(h.LegendObj(idx), 'XData');
        yy  = get(h.LegendObj(idx), 'YData');
        set(h.LegendObj(idx)  ,'Visible','off')
        set(h.LegendObj(idx+1),'Visible','off')
        for jLayer = 1 : nLayers
            if nLayers > 1
                jColor = o.LineColor(j,:) + (1 - o.LineColor(j,:)) *...
                         o.Brightness * (jLayer-1)/(nLayers-1);
            else
                jColor = o.LineColor(j,:);
            end
            jSize = 1 - (jLayer - 1)/nLayers;
            h.LegendLines(jLayer,j) = line('XData'     , xx,...
                                           'YData'     , yy,...
                                           'Parent'    , h.Legend,...
                                           'LineStyle' , o.LineStyle{j},...
                                           'Color'     , jColor,...
                                           'Marker'    , o.LineMarker{j},...
                                           'MarkerSize', o.LineMarkerSize(j)*jSize,...
                                           'LineWidth' , o.LineWidth(j)*jSize,...
                                           'MarkerFaceColor', jColor);
        end
    end
else
    h.Legend    = [];
    h.LegendObj = [];
end

%% 7. Postamble
if ~InitHold
  hold off
end

if nargout > 0
    h.Options    = o;
    varargout{1} = h;
end

%% ------------------------------------------------------------------------
