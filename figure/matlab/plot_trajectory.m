% ----------------------------------------------------------------------
% Plotting of a set of trajectory, coloring with additional field
% ----------------------------------------------------------------------

function plot = plot_trajectory (tra,varargin)

% Set plotting mode
mode=3;

% Handle the input arguments
if ( length(varargin) == 0 ) 
    col.name = 'none';
    pts      = [];
elseif ( length(varargin) == 1 )
    col      = varargin{1};
    pts      = [];
elseif ( length(varargin) == 2 )
    col      = varargin{1};
    pts      = varargin{2};
else
    error('Invalid call to plot_trajectory...');
end

% Set the selection of the trajectories
fn = fieldnames(tra);
if ( any ( strcmp(fn,'select') ) )
    tralist = tra.select;
else
    tralist = (1:tra.ntra)';
end

% Set some default parameters
fn = fieldnames(col);
if ( ~ any ( strcmp(fn,'shift') ) )
    col.shift = 47;
end
if ( ~ any ( strcmp(fn,'name') ) )
    col.name = 'col_table.txt';
end
if ( ~ any ( strcmp(fn,'orientation') ) )
    col.orientation = 'vert';
end
if ( ~ any ( strcmp(fn,'field') ) )
    col.field = 'p';
end
if ( ~ any ( strcmp(fn,'spacing') ) )
    fld = tra.(char(col.field));
    col.spacing = round(linspace(min(fld(:)),max(fld(:)),10));
end

% Get the orientation (order of values) of the colormap
n    = length(col.spacing);
diff = col.spacing(1:n-1)-col.spacing(2:n);
if ( all(diff >= 0) )
    order = 'reverse';
elseif ( all(diff <= 0) )
    order = 'normal';
else
    error('Invalid color spacing...');
end

% Set colormap
if ( ~ strcmp(col.name,'none') ) 
     c_map = col_scale(col.spacing,col.spacing);
     col_load(col.name,c_map.xtick,col.shift);
     c_map.colors = colormap;
     c_map.ncol = length(c_map.colors)
end
     
% ------- Plotting mode 1 -----------------------------------------------------------
if ( mode == 1) 
    
  for j=tralist'
 
  % Extract a single trajectory  
  tim = tra.time              (tra.label==j);
  lat = tra.lat               (tra.label==j);
  lon = tra.lon               (tra.label==j);
  fld = tra.(char(col.field)) (tra.label==j);
    
  % Draw the trajectory (either as a black line or colored with a field)
  for k=2:tra.ntime
    if ( abs(lon(k-1)-lon(k)) < 180 ) 
        
      if ( strcmp(col.name,'none') ) 
         [h]=linem([ lat(k-1) lat(k) ],[ lon(k-1) lon(k) ],'k','Linewidth',  1.0);
      else
         [h]=linem([ lat(k-1) lat(k) ],[ lon(k-1) lon(k) ],'k','Linewidth',  1.0);
         if ( strcmp(order,'normal') ) 
     	    icol=sum( 0.5*(fld(k-1)+fld(k)) > col.spacing );
         else
            icol=sum( 0.5*(fld(k-1)+fld(k)) < col.spacing );
         end
         if ( (icol>=1) & (icol <= c_map.ncol) )
	        set(h,'color',[c_map.colors(icol,1) c_map.colors(icol,2) c_map.colors(icol,3)]);
         else
            set(h,'color','k');
            set(h,'Linestyle',':');
         end
      end
      
    end
  end

  % Set time markers
  for k=1:length(pts)
     mask = (tim == pts(k) );
     if ( any(mask) ) 
        linem([ lat(mask) ],[ lon(mask) ],'marker','o', ...
              'markersize',5,'color','k','MarkerEdgeColor',[.4 .4 .4], ...
	          'MarkerFaceColor','k');
     end
  end
    
  end

end

% ------- Plotting mode 2 -----------------------------------------------------------
if ( mode == 2) 
    
  for j=tralist'
      
  % Extract a single trajectory  
  tim = tra.time              (tra.label==j);
  lat = tra.lat               (tra.label==j);
  lon = tra.lon               (tra.label==j);
  fld = tra.(char(col.field)) (tra.label==j);
    
  % Extract coordinates 
  lat1 = lat(2:tra.ntime);
  lat2 = lat(1:(tra.ntime-1));
  lon1 = lon(2:tra.ntime);
  lon2 = lon(1:(tra.ntime-1));

  % Get the color code for the line segments
  if ( strcmp(col.name,'none') ) 
    icol = zeros(size(fld1));
  else
    fld1 = fld(2:tra.ntime);
    fld2 = fld(1:(tra.ntime-1));
    ccol = kron(col.spacing,ones(size(fld1)));
    cfld = kron(0.5*(fld1+fld2),ones(size(col.spacing)));
    if ( strcmp(order,'normal') ) 
       icol = sum( cfld > ccol, 2);
    else
       icol = sum( cfld < ccol, 2);
    end
    icol( (icol<1) | (icol>c_map.ncol) ) = -1;
  end    

  % Draw colored line segments (loop through all colors)
  for i=1:c_map.ncol
     mask = (icol == i);
     if ( any(mask) ) 
  %      [h]=linem( [ lat1(mask) lat2(mask) ],[ lon1(mask) lon2(mask) ],'k','Linewidth',  1.0);
        [h]=linem( [ lat1(mask) lat2(mask) ]',[ lon1(mask) lon2(mask) ]','k','Linewidth',  1.0);
        set(h,'color',[c_map.colors(i,1) c_map.colors(i,2) c_map.colors(i,3)])
     end
  end

  % Draw line sgements outside the color range
  mask = (icol == -1);
  if ( any(mask) ) 
      [h]=linem([ lat1(mask) lat2(mask) ]',[ lon1(mask) lon2(mask) ]','k','Linewidth',  1.0);
      set(h,'color','r');
      set(h,'Linestyle','-');
  end

  % Draw black line segments (if no color is needed)
  mask = (icol == 0);
  if ( any(mask) ) 
      [h]=linem([ lat1(mask) lat2(mask) ]',[ lon1(mask) lon2(mask) ]','k','Linewidth',  1.0);
  end
  
  % Set time markers
  for k=1:length(pts)
     mask = (tim == pts(k) );
     if ( any(mask) ) 
        linem([ lat(mask) ],[ lon(mask) ],'marker','o', ...
              'markersize',5,'color','k','MarkerEdgeColor',[.4 .4 .4], ...
	          'MarkerFaceColor','k');
     end
  end

  end 
end

% ------- Plotting mode 2 -----------------------------------------------------------
if ( mode == 3 ) 
  
  % Set the total number of coordinates
  n = tra.ntime * tra.ntra;  
  
  % Extract trajectories  
  tim = tra.time;
  lat = tra.lat;
  lon = tra.lon;
  fld = tra.(char(col.field));
  
  % Mark different trajectories
  gaps = (1:tra.ntra) * tra.ntime;
  lat( gaps ) = NaN;
  lon( gaps ) = NaN;
  
  % Extract coordinates 
  lat1 = lat(2:n);
  lat2 = lat(1:(n-1));
  lon1 = lon(2:n);
  lon2 = lon(1:(n-1));

  % Get the color code for the line segments
  if ( strcmp(col.name,'none') ) 
    icol = zeros(size(fld1));
  else
    fld1 = fld(2:n);
    fld2 = fld(1:(n-1));
    ccol = kron(col.spacing,ones(size(fld1)));
    cfld = kron(0.5*(fld1+fld2),ones(size(col.spacing)));
    if ( strcmp(order,'normal') ) 
       icol = sum( cfld > ccol, 2);
    else
       icol = sum( cfld < ccol, 2);
    end
    icol( (icol<1) | (icol>c_map.ncol) ) = -1;
  end    

  % Draw colored line segments (loop through all colors)
  for i=1:c_map.ncol
     mask = (icol == i);
     if ( any(mask) ) 
        [h]=linem( [ lat1(mask) lat2(mask) ]',[ lon1(mask) lon2(mask) ]','k','Linewidth',  1.0);
        set(h,'color',[c_map.colors(i,1) c_map.colors(i,2) c_map.colors(i,3)])
     end
  end

  % Draw line sgements outside the color range
  mask = (icol == -1);
  if ( any(mask) ) 
      [h]=linem([ lat1(mask) lat2(mask) ]',[ lon1(mask) lon2(mask) ]','k','Linewidth',  1.0);
      set(h,'color','r');
      set(h,'Linestyle','-');
  end

  % Draw black line segments (if no color is needed)
  mask = (icol == 0);
  if ( any(mask) ) 
      [h]=linem([ lat1(mask) lat2(mask) ]',[ lon1(mask) lon2(mask) ]','k','Linewidth',  1.0);
  end
  
  % Set time markers
  for k=1:length(pts)
     mask = (tim == pts(k) );
     if ( any(mask) ) 
        n = length(lat(mask));
        lat3( 1:(2*n) ) = NaN; lat3( 1:2:(2*n) ) = lat(mask);
        lon3( 1:(2*n) ) = NaN; lon3( 1:2:(2*n) ) = lon(mask);
        linem([ lat3 ],[ lon3 ],'marker','o', ...
              'markersize',5,'color','k','MarkerEdgeColor',[.4 .4 .4], ...%
	          'MarkerFaceColor','k');
     end
  end
end

% ------- Plotting colorbar -----------------------------------------------------------

if ( ~ strcmp(col.name,'none') ) 
   caxis(c_map.caxis);
   if ( strcmp(col.orientation,'hori') ) 
      q=colorbar('hori');
      set(q,'xtick',c_map.xtick,'XTickLabel',c_map.label);
   elseif ( strcmp(col.orientation,'vert') ) 
      q=colorbar('vert');
      set(q,'ytick',c_map.xtick,'YTickLabel',c_map.label);
   end
end

% Return status
plot = 1;
