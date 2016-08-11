
function c_map = col_scale (dat_values,d_map);
  
% Handling of irregularly spaced color tables

% Init the label strings and output values
for i=1:length(dat_values)
  c_map.label(i) = cellstr(num2str(dat_values(i))); 
end
c_map.values=dat_values;
c_map.xtick=linspace(0,1,length(dat_values));
c_map.caxis=[ min(c_map.xtick) max(c_map.xtick) ];

% Get the transformation 
%pp=spline(dat_values,[ 0 c_map.xtick 0]);
pp = interp1(dat_values,c_map.xtick,'linear','pp');

% Scale the output array
c_map.data=d_map;
c_map.data( c_map.data < min(dat_values) ) = min(dat_values);
c_map.data( c_map.data > max(dat_values) ) = max(dat_values);
c_map.data = ppval(c_map.data,pp);
c_map.testx = linspace(min(c_map.values),max(c_map.values),100);
c_map.testy = ppval(c_map.testx,pp);




