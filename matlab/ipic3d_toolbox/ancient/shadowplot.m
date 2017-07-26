function varargout = shadowplot(varargin)
% SHADOWPLOT    Add a shadow to an existing surface or patch plot
%
% For some surface plots, it can be helpful to visualize the shadow (2D
% projection) of the surface.  This can give a quick perspective on the
% data's variance.
%
% SHADOWPLOT PLANE   Adds a shadow plot on the PLANE boundary
%  PLANE can be:
%    x, y, or z: Plots on back/top wall of x, y or z
%    1 .. 6    : Plots on Nth wall, numbered as in AXIS:
%                [xmin xmax ymin ymax zmin zmax]
%
% SHADOWPLOT(HAX,PLANE) Adds a shadow plot on the Nth wall on axes HAX
% HS = SHADOWPLOT(...) Returns a handle to the shadow (a patch)
%
% Examples:
%    figure
%    surf(peaks)
%    shading interp
%    shadowplot x   % Back X Wall
%    shadowplot y   % Back Y Wall
%
%    figure
%    surf(peaks);hold on
%    surf(peaks+10)
%    shading interp
%    hs = shadowplot(1);
%    set(hs,'FaceColor','r');     % Red shadow
%    alpha(hs,.1)                 % More transparent
%    set(hs(1),'XData',get(hs(1),'XData')*.9)  % Move farther away
%
% UPDATE (9/07): Now includes limited support for data encapsulated in
% HGTRANSFORMS, thanks to Patrick Barney (psbarne@sandia.gov).  


% Scott Hirsch
% shirsch@mathworks.com
% Copyright 2004-2007 The MathWorks, Inc

%% We define three dimensions.  1=x, 2=y, 3=z
% dimplane - dimension that's constant in the projection plane (user-specified)
% dimvar - dimension in which data varies (typically 3)
% dimother - the other dimension (couldn't come up with a good name!).  

%% Parse input arguments.
if nargin==1
    hAx = gca;
    plane = lower(varargin{1});
elseif nargin==2
    hAx = varargin{1};
    plane = lower(varargin{2});
end;

%% Convert plane to numeric dimension
% plane can be specified as a string (x,y,z) or as a number (1..6)
if ~isstr(plane)
    dimplane = ceil(plane/2);
    axind = plane;  % Index into AXIS to get boundary plane
else    % string
    switch plane
        case 'x'
            dimplane = 1;
            axind = 2;  % Index into AXIS to get boundary plane
        case 'y'
            dimplane = 2;
            axind = 4;
        case 'z'
            dimplane = 3;
            axind = 6;
        otherwise
            error('Plane must be one of: ''x'', ''y'', or ''z'' or a number between 1 and 6');
    end;
end;

%% Get coordinates for placing surface from axis limits
ax = axis;
% ============ force axis into 3d mode =============
if length(axis==4)
  % axis problem. get the current view, rotate it, then
  % redo the axis and return to the original view.
  [az,el] = view;
  view(45,45)
  ax = axis;
  view(az,el)
end
planecoord = ax(axind);     % Plane Coordinate - back wall

%% Turn hold on
hold_current = ishold(hAx);
if hold_current == 0
    hold_current = 'off';
else
    hold_current = 'on';
end;

hold(hAx,'on')

%% Get handles to all surfaces
kids = findobj(hAx,'Type','surface');
h = [];

% Also get handles to all patch objects
kidsp = findobj(hAx,'Type','patch');
hp = [];

for ii=1:length(kids)       % Do separately for each surface
    hSurf = kids(ii);   % Current surface

    % We do everything with the X, Y, and ZData of the surface
    surfdata = get(hSurf,{'XData','YData','ZData'});

    % XData and YData might be vectors or matrices.  Force them to be
    % matrices (a la griddata)
    [Ny,Nx] = size(surfdata{3});
    if isvector(surfdata{1})
        surfdata{1} = repmat(surfdata{1},Ny,1);
    end;
    if isvector(surfdata{2})
        surfdata{2} = repmat(surfdata{2},1,Nx);
    end;

    % Figure out which two axes are independent (i.e., monotonic)
    grids = [ismeshgrid(surfdata{1}) ismeshgrid(surfdata{2}) ismeshgrid(surfdata{3})];
    if sum(grids)<2, error('Surface must have at least 2 monotonically increasing dimensions');end

    % The remaining dimension is the one along which data varies
    dimvar  = find(~grids); % Dimension where data varies
    if isempty(dimvar)  % All 3 dimensions are monotonic.  not sure what to do
        dimvar = max(setdiff(1:3,dimplane));% pick largest value that isn't dimplane
    end;
    
    if dimvar==dimplane
        error('Can not project data in the dimension that varies.  Try another plane')
    end;

    %dimdiff: dimension for taking difference (figure out through trial and error)
    % dimplane=1, dimvar=3: 2
    % dimplane=1, dimvar=2: 2
    % dimplane=2, dimvar=1: 2
    % dimplane=2, dimvar=3: 1
    % dimplane=3, dimvar=2: 1
    % dimplane=3, dimvar=1: 1
    dimdiff = 2;    % Most cases
    if (dimplane==2&&dimvar==3) | (dimplane==3)
        dimdiff = 1;
    end;
    
    % Compute projection
    dmin = min(surfdata{dimvar},[],dimdiff);    % Min of data projected onto this plane
    dmax = max(surfdata{dimvar},[],dimdiff);    % Max of data projected onto this plane

    dmin = dmin(:); % Force into row vector
    dmax = dmax(:);

    nval = length(dmin)*2 + 1;  % Total number of values we'll use for shadow

    % Compute shadow coordinates
    % Pull out independent variable
    dimother = setxor([dimvar dimplane],1:3); % Remaining dimension
    d1 = surfdata{dimother}(:,1);    % Not sure if should take row or col. find the dimension that varies
    d2 = surfdata{dimother}(1,:);
    if d1(1) ~= d1(end)
        dind = d1;
    else
        dind = d2';
    end;


    shadow{dimplane} = repmat(planecoord,nval,1);       % In the plane
    shadow{dimother} = [dind;flipud(dind);dind(1)];     % Independent variable
    shadow{dimvar} = [dmin;flipud(dmax);dmin(1)];       % the varying data

    h(ii) = patch(shadow{1},shadow{2},shadow{3},[.3 .3 .3]);
    alpha(h(ii),.3)
    set(h(ii),'LineStyle','none')
    % set a tag, so that a shadow will not try to cast a shadow
    set(h(ii),'Tag','Shadow')
end;

%% Shadow any patches, unless they are already shadows.
hp = [];

for ii=1:length(kidsp)       % Do separately for each patch object
    hPat = kidsp(ii);   % Current patch

    % Is this patch already tagged as a Shadow?
    if ~strcmpi(get(hPat,'Tag'),'Shadow')

      % We do everything with the X, Y, and ZData of the surface
      patdata = get(hPat,{'XData','YData','ZData'});

      switch get(get(hPat,'par'),'Type')
          case 'hgtransform'
              M=get(get(hPat,'par'),'Matrix');
              try
                  switch get(get(get(hPat,'par'),'par'),'type')
                      case 'hgtransform'
                          M2=get(get(get(hPat,'par'),'par'),'Matrix');
                          M=M2*M;
                  end
              end
              M=M(1:3,1:3);
              xyz=[patdata{1}(:),patdata{2}(:),patdata{3}(:)]*M';
              [n,m]=size(patdata{1});
              patdata{1}=reshape(xyz(:,1),n,m);
              patdata{2}=reshape(xyz(:,2),n,m);
              patdata{3}=reshape(xyz(:,3),n,m);
          otherwise
      end
      % Just replace the x, y, or z coordinate as indicated by dimplane
      patdata{dimplane} = repmat(planecoord,size(patdata{dimplane}));

      % then its just a call to patch
      hp(ii) = patch(patdata{1},patdata{2},patdata{3},[.3 .3 .3]);
      alpha(hp(ii),.3)
      set(hp(ii),'LineStyle','none')
      
      % set a tag, so that a shadow will not try to cast a shadow
      set(hp(ii),'Tag','Shadow')
    end
end;
h=[h,hp];

hold(hAx,hold_current)  % Return to original state

if nargout
    varargout{1} = h;
end;

function isgrid = ismeshgrid(d)
% Check if d looks like it came from griddata
dd = diff(d);
ddt = diff(d');

if all(~dd(:)) | all(~ddt(:))
    isgrid = 1;
else
    isgrid = 0;
end;
% if ~any(d(:,1) - d(:,end)) | ~any(d(1,:) - d(end,:))
%     isgrid = 1;
% else
%     isgrid = 0;
% end;


