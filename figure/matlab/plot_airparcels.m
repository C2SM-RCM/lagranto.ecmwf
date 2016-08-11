% ------------------------------------------------------------------------
% Preparations
% ------------------------------------------------------------------------

% Add some matlab paths
addpath('/usr/local/matlabtools/mexnc/');
addpath('/usr/local/matlabtools/snctools/');
addpath('/usr/local/matlabtools/cdf_io/');
addpath('/home/sprenger/lagranto/matlab/');

% Set the datadirectories
cdfdir   = '/lhome/sprenger/lagranto/test/cdf/';
tradir   = '/home/sprenger/lagranto/matlab/';

% Read a trajectory file
tra = read_trajectory( [ tradir '/trajectory2.lsl' ]);

% Set the tabel for input files
filelist.n = 4;
filelist.date(1) = cellstr('19891020_00'); filelist.time(1) =  0;
filelist.date(2) = cellstr('19891020_06'); filelist.time(2) =  6;
filelist.date(3) = cellstr('19891020_12'); filelist.time(3) = 12;
filelist.date(4) = cellstr('19891020_18'); filelist.time(4) = 18; 

% Get the grid parameters from first data file
inp                     = ncget([ cdfdir 'P' char(filelist.date(1)) ],'U','PS');
grid.nx                 = size(inp.U.data,3);
grid.ny                 = size(inp.U.data,2);
grid.nz                 = size(inp.U.data,1);
grid.xmin               = inp.cstfile.lonmin;
grid.ymin               = inp.cstfile.latmin;
grid.dx                 = inp.cstfile.dellon;
grid.dy                 = inp.cstfile.dellat;
grid.xmax               = grid.xmin + double(grid.nx-1) * grid.dx;
grid.ymax               = grid.ymin + double(grid.ny-1) * grid.dy;
grid.mdv                = -999.
grid.lon1               = grid.xmin + double(0:grid.nx-1) * grid.dx;
grid.lat1               = grid.ymin + double(0:grid.ny-1) * grid.dy;
[ grid.lon2 grid.lat2 ] = meshgrid(grid.lon1,grid.lat1);
grid.aklay              = inp.cstfile.aklay;
grid.aklev              = inp.cstfile.aklev;
grid.bklay              = inp.cstfile.bklay;
grid.bklev              = inp.cstfile.bklev;

% ------------------------------------------------------------------------
% Loop over al trajectory times
% ------------------------------------------------------------------------

% Set flags for load manager
load0 = -1;
load1 = -1;

% No trajectories are selected
select = zeros(tra.ntra,1);

% Init the textlines
for j=1:tra.nfield
    textline(j) = cellstr(' ');
end
first = 1;

% Loop over all trajectory times
i     = 1;
stat  = 1;
while ( stat ~= 0 )

    % Get the time 
    time    = tra.time(i);
    
    % Decide which data files are needed
    need0 = find( time <= filelist.time , 1 );
    need1 = need0 + 1;
    if ( need1 > filelist.n ) 
       need1 = need0;
    end
       
    % Load files if not yet done
    if ( need0 == load1 )
        inp0  = inp1;
        load0 = need0;
    end
    if ( need1 == load0 )
        inp1  = inp0;
        load1 = need1;
    end
    if ( need0 ~= load0 ) 
        inps        = ncget([ cdfdir 'S' char(filelist.date(need0)) ],'TH','PV','RH');
        inpp        = ncget([ cdfdir 'P' char(filelist.date(need0)) ],'Q','U','V','T','PS');
        inp0        = inpp;
        inp0.TH     = inps.TH;
        inp0.PV     = inps.PV;
        inp0.RH     = inps.RH;
        inp0.P      = inp.PS;
        inp0.P.data = kron(inp0.PS.data,grid.bklay) + kron(ones(grid.ny,grid.nx),grid.aklay);
        inp0.P.data = reshape(inp0.P.data,grid.nz,grid.ny,grid.nx);
        load0       = need0;
    end
    if ( need1 ~= load1 ) 
        inps        = ncget([ cdfdir 'S' char(filelist.date(need1)) ],'TH','PV','RH');
        inpp        = ncget([ cdfdir 'P' char(filelist.date(need1)) ],'Q','U','V','T','PS');
        inp1        = inpp;
        inp1.TH     = inps.TH;
        inp1.PV     = inps.PV;
        inp1.RH     = inps.RH;
        inp1.P      = inp.PS;
        inp1.P.data = kron(inp1.PS.data,grid.bklay) + kron(ones(grid.ny,grid.nx),grid.aklay);
        inp1.P.data = reshape(inp1.P.data,grid.nz,grid.ny,grid.nx);
        load1       = need1;
    end

    % Get the mean pressure of the trajectory sample
    if ( any(select) )
        mask = ( find(select) - 1) * tra.ntime + i;
    else
        mask = ( (1:tra.ntra) -1 ) * tra.ntime + i;
    end
    hori.level = mean ( tra.p ( mask ) );
    
    % Interpolate fields to a pressure surface
    [ xpos  ypos zpos ] = meshgrid( grid.lon1, grid.lat1, hori.level );
    ind                 = get_index(inp0.P.data,[ zpos(:) ypos(:) xpos(:) ],grid);
    hori0.U             = reshape( int_index(inp0.U.data,ind,grid.mdv), 1,grid.ny,grid.nx);
    hori0.V             = reshape( int_index(inp0.V.data,ind,grid.mdv), 1,grid.ny,grid.nx);
    hori0.VEL           = sqrt(hori0.U .^ 2 + hori0.V .^ 2);

    [ xpos  ypos zpos ] = meshgrid( grid.lon1, grid.lat1, hori.level );
    ind                 = get_index(inp1.P.data,[ zpos(:) ypos(:) xpos(:) ],grid);
    hori1.U             = reshape( int_index(inp1.U.data,ind,grid.mdv), 1,grid.ny,grid.nx);
    hori1.V             = reshape( int_index(inp1.V.data,ind,grid.mdv), 1,grid.ny,grid.nx);
    hori1.VEL           = sqrt(hori1.U .^ 2 + hori1.V .^ 2);

    % Do a time interpolation (and get rid of 1-dimensions) 
    if ( need1 ~= need0 )
        frac     = ( time - filelist.time(need0) ) / ( filelist.time(need1) - filelist.time(need0) );
        hori.U   = squeeze( (1-frac) * hori0.U   + frac * hori1.U   );
        hori.V   = squeeze( (1-frac) * hori0.V   + frac * hori1.V   );
        hori.VEL = squeeze( (1-frac) * hori0.VEL + frac * hori1.VEL );
    else
        hori.U   = squeeze( hori0.U   );
        hori.V   = squeeze( hori0.V   );
        hori.VEL = squeeze( hori0.VEL );
    end
    
    % Open a new figure ans set the geographical projection
    figure(1);
    clf;
    load coast
    %h=axesm('MapProjection','stereo','origin',[ 90 0 ]);
    h=axesm( 'eqdcylin','MapLatLimit',[ 20 80 ], 'MapLonLimit', [ -100 40 ] );
    gridm;
    h=plotm(lat,long,'Color','k','LineWidth',1.5);
    %axis([-0.6 0.67 -1.2 -0.1]);
    axis([ -1.0    1.0    0.4000    1.4000 ]);
    daspect([ 1 0.75 1] );
    title ( [ 'P = ' num2str(round(hori.level)) ' hPa / T = ' num2str(time) ' h'],...
           'FontSize',22)

    % Plot velocity and overlay wind vectors
    c_map = col_scale([0:2.5:50 ],hori.VEL);
    [C,h] = contourfm(grid.lat1,grid.lon1,c_map.data,c_map.xtick);
    colormap('default');
    col_load( [ pwd '/col_table.txt' ],c_map.xtick,44); 
    caxis(c_map.caxis)
    q=colorbar('vert');
    set(q,'ytick',c_map.xtick,'YTickLabel',c_map.label,'FontSize',15); 
    [ c_map.x c_map.y ] = meshgrid(grid.lon1,grid.lat1);
    c_map.u             = squeeze(hori.U);
    c_map.v             = squeeze(hori.V);
    maskx               = 1:2:size(grid.lon2,2);        
    masky               = 1:2:size(grid.lon2,1);
    hold on
    q=quiverm(grid.lat2(masky,maskx),grid.lon2(masky,maskx),...
              hori.V   (masky,maskx),hori.U   (masky,maskx),'k');
    set(q,'Linewidth',1,'Color','k')
    hold off  
    
    % Handle user commands
    stat = 2;
    while ( stat == 2 )
        
        % Show the positions of the trajectory points
        hold on
        for j=1:tra.ntra
            if select(j) 
                mask = ( tra.time == time ) & ( tra.label == j );
                plotm(tra.lat(mask),tra.lon(mask),'o','markersize',10,'color','r',...
                      'MarkerEdgeColor','k','MarkerFaceColor','r');
            else
                mask = ( tra.time == time ) & ( tra.label == j );
                plotm(tra.lat(mask),tra.lon(mask),'o','markersize',10,'color','k',...
                      'MarkerEdgeColor','k','MarkerFaceColor','k');
            end
        end
        hold off
      
        % Get the cursor position
        [ x y ] = ginput(1)
        plotdom = axis;

       % Step backward
       if ( x < plotdom(1) ) 
           i = i - 1;
           if ( i < 1)
               i = tra.ntime;
           end
           stat = 1;
           
       % Step forward
       elseif ( x > plotdom(2) ) 
           i = i + 1;
           if ( i > tra.ntime )
               i = 1;
           end
           stat = 1;
           
       % Exit
       elseif ( y > plotdom(4) )
           stat = 0;
        
       % Write trajectory info
       elseif ( (y < plotdom(3)) & any(select) )
           for j=1:tra.nfield
              mask = ( find(select) - 1) * tra.ntime + i;
              val  = tra.(char(tra.field(j)));
              val  = mean( val(mask) );
              if ( first == 0 ) 
                 text(plotdom(1),plotdom(3)-0.04*j,char(textline(j)),...
                      'FontSize',15,'Color','g' );
              end
              textline(j) = cellstr([ char(tra.field(j)) ' = ' num2str(val)]);
              text(plotdom(1),plotdom(3)-0.04*j,char(textline(j)),...
                   'FontSize',15,'Color','b' );
              first = 0;
           end
        
       % Select trajectory points
       else
           [ lat0, lon0 ] = minvtran(x,y);
           mask           = (tra.time == time );
           lon1           = tra.lon(mask);
           lat1           = tra.lat(mask);
           dist           = distance(lat0,lon0,lat1,lon1);
           label          = find( min(dist) == dist );
           if  ( dist(label) < 5 ) 
              if ( select(label) == 0 ) 
                 select(label) = 1;
              elseif ( select(label) == 1 )
                 select(label) = 0;
              end
           else
               stat = 1;
           end
       end
       
    end

% End loop 
end

