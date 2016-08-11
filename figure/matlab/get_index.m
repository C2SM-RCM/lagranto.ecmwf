function out = get_index (vert,position,grid);

% General index calculation for data grid
% Usage:   out = int_index(vert,coords,grid)   
% Example: grid.xmin = -180;
%          grid.dx   =    1;
%          grid.nx   =  361;
%          grid.ymin =    0;
%          grid.dy   =    1;
%          grid.ny   =   91;
%          out get_index(inp.P.data,[ 500 45 78; 250 56.7 78.4],grid)   

% Get the position arrays
xpos = position(:,3)';
ypos = position(:,2)';
zpos = position(:,1)';

% Get the grid dimension 
nx = size(vert,3);
ny = size(vert,2);
nz = size(vert,1);

% Get the horiozontal grid indices
rindx = (xpos-grid.xmin)/grid.dx+1;
rindy = (ypos-grid.ymin)/grid.dy+1;

% Special treatment if the the <vert> field is two-dimensional
if ( nz == 1 )
   out = [ ones(size(rindy))', rindy', rindx' ];
   return;
end 
    
% Direction of the vertical axis (Pressure: descencding, Height: ascending)
dir = vert(2:nz,:,:) - vert(1:nz-1,:,:);
if ( mean(dir(:)) > 0 ) 
     zpos = -zpos;
     vert = -vert;
end

% Get the weights for the individual grid points
fracx = 1-rindx+floor(rindx);
fracy = 1-rindy+floor(rindy);

% Combine the indices i and j to a common index ij
vert1 = reshape(vert,nz,nx*ny);

% The position array and the vertical index is needed in the same format
zpos1 = kron(ones(nz,1),zpos);
zind1 = kron([1:nz]',ones(size(zpos)));
shif1 = (0:length(zpos)-1) * nz;

% Get the vertical indices at position 00
indx           = floor(rindx);
indy           = floor(rindy);
indxy          = indy+(indx-1)*ny;
isnok1         = (indx<1) | (indx>nx) | (indy<1) | (indy>ny);
indxy(isnok1)  = 1;  
vert           = vert1(:,indxy);
indz           = max((vert>=zpos1).*zind1);
isnok2         = (indz<1) | (indz>=nz);
indz(isnok2)   = 1;
indu           = indz     + shif1;
indo           = indz + 1 + shif1;
diff           = vert(indo)-vert(indu);
isnok3         = (diff==0);
diff(isnok3)   = 1;
rind00         = indz + (zpos-vert(indu)) ./ diff;
isnok4         = isnok1 | isnok2 | isnok3;
rind00(isnok4) = NaN;

% Get the vertical indices at positions 01, 10, 11 (if necessary)
if ( any(fracx ~= 1) | any(fracy ~=1) ) 
    
    indx           = floor(rindx+1);
    indy           = floor(rindy);
    indxy          = indy+(indx-1)*ny;
    isnok1         = (indx<1) | (indx>nx) | (indy<1) | (indy>ny);
    indxy(isnok1)  = 1;  
    vert           = vert1(:,indxy);
    indz           = max((vert>zpos1).*zind1);
    isnok2         = (indz<1) | (indz>=nz);
    indz(isnok2)   = 1;
    indu           = indz     + shif1;  
    indo           = indz + 1 + shif1;
    diff           = vert(indo)-vert(indu);
    isnok3         = (diff==0);
    diff(isnok3)   = 1;
    rind01         = indz + (zpos-vert(indu)) ./ diff;
    isnok4         = isnok1 | isnok2 | isnok3;
    rind01(isnok4) = NaN;

    indx           = floor(rindx);
    indy           = floor(rindy+1);
    indxy          = indy+(indx-1)*ny;
    isnok1         = (indx<1) | (indx>nx) | (indy<1) | (indy>ny);
    indxy(isnok1)  = 1;  
    vert           = vert1(:,indxy);
    indz           = max((vert>zpos1).*zind1);
    isnok2         = (indz<1) | (indz>=nz);
    indz(isnok2)   = 1;
    indu           = indz     + shif1;  
    indo           = indz + 1 + shif1;
    diff           = vert(indo)-vert(indu);
    isnok3         = (diff==0);
    diff(isnok3)   = 1;
    rind10         = indz + (zpos-vert(indu)) ./ diff;
    isnok4         = isnok1 | isnok2 | isnok3;
    rind10(isnok4) = NaN;

    indx           = floor(rindx+1);
    indy           = floor(rindy+1);
    indxy          = indy+(indx-1)*ny;
    isnok1         = (indx<1) | (indx>nx) | (indy<1) | (indy>ny);
    indxy(isnok1)  = 1;  
    vert           = vert1(:,indxy);
    indz           = max((vert>zpos1).*zind1);
    isnok2         = (indz<1) | (indz>=nz);
    indz(isnok2)   = 1;
    indu           = indz     + shif1;  
    indo           = indz + 1 + shif1;
    diff           = vert(indo)-vert(indu);
    isnok3         = (diff==0);
    diff(isnok3)   = 1;
    rind11         = indz + (zpos-vert(indu)) ./ diff;
    isnok4         = isnok1 | isnok2 | isnok3;
    rind11(isnok4) = NaN;

else

    rind01 = rind00;
    rind10 = rind00;
    rind11 = rind00;
    
end
    
% Vectorised interpolation
rindz = fracx     .* fracy     .* rind00 + ...
        (1-fracx) .* fracy     .* rind01 + ...
        fracx     .* (1-fracy) .* rind10 + ...
        (1-fracx) .* (1-fracy) .* rind11;
        
% Set the output of the function
nok = isnan(rindz') | isnan(rindy') | isnan(rindx');
out = [ rindz', rindy', rindx' ];



