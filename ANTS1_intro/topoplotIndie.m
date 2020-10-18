function [handle,pltchans,epos,Zi] = topoplotIndie(Values,chanlocs,varargin)

%% Set defaults

headrad = 0.5;          % actual head radius - Don't change this!
GRID_SCALE = 67;        % plot map on a 67X67 grid
CIRCGRID   = 201;       % number of angles to use in drawing circles
HEADCOLOR = [0 0 0];    % default head color (black)
HLINEWIDTH = 1.7;         % default linewidth for head, nose, ears
BLANKINGRINGWIDTH = .035;% width of the blanking ring
HEADRINGWIDTH    = .007;% width of the cartoon head ring
plotrad = .6;
Values = double(Values);
SHADING = 'interp';
CONTOURNUM = 6;
ELECTRODES = 'on';

nargs = nargin;
if nargs > 2
    for i = 1:2:length(varargin)
        Param = lower(varargin{i});
        Value = varargin{i+1};
        switch Param
            case 'numcontour'
                CONTOURNUM = Value;
            case 'electrodes'
                ELECTRODES = lower(Value);
            case 'plotrad'
                plotrad = Value;
            case 'shading'
                SHADING = lower(Value);
                if ~any(strcmp(SHADING,{'flat','interp'}))
                    error('Invalid shading parameter')
                end
        end
    end
end

Values = Values(:); % make Values a column vector

%% Read channel location

labels = {chanlocs.labels};
Th = [chanlocs.theta];
Rd = [chanlocs.radius];

Th = pi/180*Th;                              % convert degrees to radians
allchansind = 1:length(Th);
plotchans = 1:length(chanlocs);
[x,y] = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates

%% remove infinite and NaN values

inds = union(find(isnan(Values)), find(isinf(Values))); % NaN and Inf values
for chani=1:length(chanlocs)
    if isempty(chanlocs(chani).X); inds = [inds chani]; end
end

plotchans   = setdiff(plotchans,inds);

plotchans   = abs(plotchans);   % reverse indicated channel polarities
allchansind = allchansind(plotchans);
Th          = Th(plotchans);
Rd          = Rd(plotchans);
x           = x(plotchans);
y           = y(plotchans);
labels      = char(labels(plotchans)); % remove labels for electrodes without locations
Values      = Values(plotchans);
intrad      = min(1.0,max(Rd)*1.02);             % default: just outside the outermost electrode location

%% Find plotting channels

pltchans = find(Rd <= plotrad); % plot channels inside plotting circle
intchans = find(x <= intrad & y <= intrad); % interpolate and plot channels inside interpolation square

%% Eliminate channels not plotted

allx  = x;
ally  = y;
allchansind = allchansind(pltchans);
intTh = Th(intchans);           % eliminate channels outside the interpolation area
intRd = Rd(intchans);
intx  = x(intchans);
inty  = y(intchans);
Th    = Th(pltchans);              % eliminate channels outside the plotting area
Rd    = Rd(pltchans);
x     = x(pltchans);
y     = y(pltchans);

intValues = Values(intchans);
Values = Values(pltchans);

labels= labels(pltchans,:);

%% Squeeze channel locations to <= headrad
squeezefac = headrad/plotrad;
intRd = intRd*squeezefac; % squeeze electrode arc_lengths towards the vertex
Rd    = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex
% to plot all inside the head cartoon
intx  = intx*squeezefac;
inty  = inty*squeezefac;
x     = x*squeezefac;
y     = y*squeezefac;
allx  = allx*squeezefac;
ally  = ally*squeezefac;

%% create grid
xmin = min(-headrad,min(intx)); xmax = max(headrad,max(intx));
ymin = min(-headrad,min(inty)); ymax = max(headrad,max(inty));
xi   = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
yi   = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)

[Xi,Yi,Zi] = griddata(inty,intx,intValues,yi',xi,'v4'); % interpolate data

%% Mask out data outside the head

mask = (sqrt(Xi.^2 + Yi.^2) <= headrad); % mask outside the plotting circle
Zi(mask == 0) = NaN;                  % mask non-plotting voxels with NaNs
grid = plotrad;                       % unless 'noplot', then 3rd output arg is plotrad
delta = xi(2)-xi(1); % length of grid entry

%% Scale the axes and make the plot
cla  % clear current axis
hold on
h = gca; % uses current axes
AXHEADFAC = 1.05;     % do not leave room for external ears if head cartoon
set(gca,'Xlim',[-headrad headrad]*AXHEADFAC,'Ylim',[-headrad headrad]*AXHEADFAC);
unsh = (GRID_SCALE+1)/GRID_SCALE; % un-shrink the effects of 'interp' SHADING

if strcmp(SHADING,'interp')
    handle = surface(Xi*unsh,Yi*unsh,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING);
else
    handle = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor',SHADING);
end
contour(Xi,Yi,Zi,CONTOURNUM,'k','hittest','off');

%% Plot filled ring to mask jagged grid boundary
hwidth = HEADRINGWIDTH;                   % width of head ring
hin  = squeezefac*headrad*(1- hwidth/2);  % inner head ring radius

if strcmp(SHADING,'interp')
    rwidth = BLANKINGRINGWIDTH*1.3;             % width of blanking outer ring
else
    rwidth = BLANKINGRINGWIDTH;         % width of blanking outer ring
end
rin    =  headrad*(1-rwidth/2);              % inner ring radius
if hin>rin
    rin = hin;                              % dont blank inside the head ring
end

circ = linspace(0,2*pi,CIRCGRID);
rx = sin(circ);
ry = cos(circ);
ringx = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
ringy = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];
ringh = patch(ringx,ringy,0.01*ones(size(ringx)),get(gcf,'color'),'edgecolor','none','hittest','off'); hold on

%% Plot cartoon head, ears, nose
headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
ringh = patch(headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR,'hittest','off'); hold on

% Plot ears and nose
base  = headrad-.0046;
basex = 0.18*headrad;                   % nose width
tip   = 1.15*headrad;
tiphw = .04*headrad;                    % nose tip half width
tipr  = .01*headrad;                    % nose tip rounding
q = .04; % ear lengthening
EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % headrad = 0.5
EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
sf    = headrad/plotrad;                                          % squeeze the model ears and nose
% by this factor
plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,2*ones(size([basex;tiphw;0;-tiphw;-basex])),'Color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off');                 % plot nose
plot3(EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')    % plot left ear
plot3(-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH,'hittest','off')   % plot right ear

%% Mark electrode locations
if strcmp(ELECTRODES,'on')   % plot electrodes as spots
    hp2 = plot3(y,x,ones(size(x)),'.','Color',[0 0 0],'markersize',5,'linewidth',.5,'hittest','off');
elseif strcmp(ELECTRODES,'labels')  % print electrode names (labels)
    for i = 1:size(labels,1)
        text(double(y(i)),double(x(i)),1,labels(i,:),'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'hittest','off')
    end
elseif strcmp(ELECTRODES,'numbers')
    for i = 1:size(labels,1)
        text(double(y(i)),double(x(i)),1,int2str(allchansind(i)),'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'hittest','off')
    end
end

epos=[x; y];
axis off
axis equal

