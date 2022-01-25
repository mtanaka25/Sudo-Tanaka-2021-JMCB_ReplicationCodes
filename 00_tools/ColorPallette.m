function ColorList = ColorPallette(varargin)
% 

%% ------------------------------------------------------------------------
%% 1. Default Options
ColorOptions.MyPallette = [... 
  [ 192,  80,  77];... % red
  [  79, 129, 189];... % blue
  [ 155, 187,  89];... % green
  [ 247, 150,  70];... % orange
  [   0,   0,   0];... % black
  [ 128, 100, 162] ... % purple
   ]./255;

ColorOptions.Matlab = [...
  [   0,   0,   1];... % blue
  [   0,  .5,   0];... % green
  [   1,   0,   0];... % red
  [ .25, .25, .25];... % grey
  [   0, .55, .55];... % dark blue
  [  .9,  .6, .05] ... % orange
   ];

ColorOptions.Excel2010 = [... 
  [ 192,  80,  77];... % red
  [  79, 129, 189];... % blue
  [ 155, 187,  89];... % green
  [ 247, 150,  70];... % orange
  [ 128, 100, 162];... % purple
  [   0,   0,   0] ... % black
   ]./255;

BaseFactor     = 1;
LightFactors   = [0,0.25,0.5];
isShowPallette = 0;
isNewFigure    = 0;
isLinePlot     = 0;
LineWidth      = 2;
UseColors      = 'MyPallette';

%% 2. Check inputs and Overwrite options
for jO = 1:(length(varargin)/2) 
  eval(sprintf('%s = varargin{%.0f};', varargin{(jO-1)*2+1}, jO*2))
end

%% 3. Make colors
BaseColors    = BaseFactor * ColorOptions.(UseColors);
nBaseColors   = size(BaseColors,1);
nLightFactors = length(LightFactors);

if ~exist('nColors','var')
    nColors = nBaseColors * nLightFactors;
end

for j = 1 : nLightFactors
    ColorList((j-1)*nBaseColors+(1:nBaseColors),:) = BaseColors + LightFactors(j) * (1 - BaseColors);
end
ColorList = ColorList(1:nColors,:);

%% 4. Plot Sample
if isShowPallette && ~isLinePlot
    if isNewFigure,figure,end
    plot(0:1, 0:1, 'w')
    hold on
    for j = 1 : nColors
        fill([j-1, j-1, j, j]/nColors,...
             [0, 0.5+j/nColors/2, 0.5+j/nColors/2, 0],...
             ColorList(j, :))
    end
    hold off
elseif  isShowPallette && isLinePlot
    if isNewFigure,figure,end
    plot(0:1, 0:1, 'w')
    hold on
    for j = 1 : nColors
        plot(0:0.1:1, j/(nColors+1)/3+(0:.1:1)*j/(nColors+1), ...
            'Color'    ,ColorList(j,:),...
            'LineWidth',LineWidth)
    end
    hold off
end


%% ------------------------------------------------------------------------
