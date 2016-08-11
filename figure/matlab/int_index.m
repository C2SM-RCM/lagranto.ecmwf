function out = int_index (field,indices,mdv)

% General interpolation routine for data grid
% Usage:   out = int_index(in,index,mdv)   
% Example: out int_index(inp.TH.data,[10.1 20.4 30.3; 23.6 56.7 78.4],mdv)   

% Check whether input data are correct
if ( length(size(field)) > 3 )
    error('Error in int_index: input field can be at most three-dimensional');
end
if ( length(size(indices)) ~= 2 | ...
     length(size(indices)) == 2 & size(indices,2) ~= length(size(field)) )
    error('Input fields for int_index are inconsistent');
end

% Get the dimensions of the input field and the indices
if ( length(size(field)) == 3 )
    nx   = size(field,3); 
    ny   = size(field,2);
    nz   = size(field,1);
    xpos = indices(:,3);
    ypos = indices(:,2);
    zpos = indices(:,1);
elseif ( length(size(field)) == 2 )
    nx   = size(field,2); 
    ny   = size(field,1);
    nz   = 1;
    xpos = indices(:,2);
    ypos = indices(:,1);
    zpos = ones(size(xpos));
end

% Make the field one-dimensional (needed for vectorised computation)
field1 = field(:);

% Get the indices of the neighbouring data points (surrounding grid box)
ix000 = max(floor(xpos),1);
iy000 = max(floor(ypos),1);
iz000 = max(floor(zpos),1);

ix001 = min(1+floor(xpos),nx);
iy001 = iy000;
iz001 = iz000;

ix010 = ix000;
iy010 = min(1+floor(ypos),ny);
iz010 = iz000;

ix011 = min(1+floor(xpos),nx);
iy011 = min(1+floor(ypos),ny);
iz011 = iz000;

ix100 = ix000;
iy100 = iy000;
iz100 = min(1+floor(zpos),nz);

ix101 = min(1+floor(xpos),nx);
iy101 = iy000;
iz101 = min(1+floor(zpos),nz);

ix110 = ix000;
iy110 = min(1+floor(ypos),ny);
iz110 = min(1+floor(zpos),nz);

ix111 = min(1+floor(xpos),nx);
iy111 = min(1+floor(ypos),ny);
iz111 = min(1+floor(zpos),nz);

% Get the weights for the individual grid points
fracx = (1-xpos+ix000);
fracy = (1-ypos+iy000);
fracz = (1-zpos+iz000);

frac000 =    fracx  .*    fracy  .*    fracz ;
frac001 = (1-fracx) .*    fracy  .*    fracz ;
frac010 =    fracx  .* (1-fracy) .*    fracz ;
frac011 = (1-fracx) .* (1-fracy) .*    fracz ;
frac100 =    fracx  .*    fracy  .* (1-fracz);
frac101 = (1-fracx) .*    fracy  .* (1-fracz);
frac110 =    fracx  .* (1-fracy) .* (1-fracz);
frac111 = (1-fracx) .* (1-fracy) .* (1-fracz);

% Vectorised interpolation 
out   = frac000 .* field1( iz000 + (iy000-1) * nz + (ix000-1) * nz * ny ) + ...
        frac001 .* field1( iz001 + (iy001-1) * nz + (ix001-1) * nz * ny ) + ...
        frac010 .* field1( iz010 + (iy010-1) * nz + (ix010-1) * nz * ny ) + ...
        frac011 .* field1( iz011 + (iy011-1) * nz + (ix011-1) * nz * ny ) + ...
        frac100 .* field1( iz100 + (iy100-1) * nz + (ix100-1) * nz * ny ) + ...
        frac101 .* field1( iz101 + (iy101-1) * nz + (ix101-1) * nz * ny ) + ...
        frac110 .* field1( iz110 + (iy110-1) * nz + (ix110-1) * nz * ny ) + ...
        frac111 .* field1( iz111 + (iy111-1) * nz + (ix111-1) * nz * ny );
    
% Handling of missing data values
mask = ( field1( iz000 + (iy000-1) * nz + (ix000-1) * nz * ny ) == mdv ) | ...
       ( field1( iz001 + (iy001-1) * nz + (ix001-1) * nz * ny ) == mdv ) | ...
       ( field1( iz010 + (iy010-1) * nz + (ix010-1) * nz * ny ) == mdv ) | ...
       ( field1( iz011 + (iy011-1) * nz + (ix011-1) * nz * ny ) == mdv ) | ...
       ( field1( iz100 + (iy100-1) * nz + (ix100-1) * nz * ny ) == mdv ) | ...
       ( field1( iz101 + (iy101-1) * nz + (ix101-1) * nz * ny ) == mdv ) | ...
       ( field1( iz110 + (iy110-1) * nz + (ix110-1) * nz * ny ) == mdv ) | ...
       ( field1( iz111 + (iy111-1) * nz + (ix111-1) * nz * ny ) == mdv );
out(mask) = mdv;    
    
    
  

    
    