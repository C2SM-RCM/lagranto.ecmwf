% ----------------------------------------------------------------------
% Plot a set of trajectories
% ----------------------------------------------------------------------

% First image (otherwise first image is not correctly written)
filename = 'test.eps';
close;
fh=figure('Units','pixels','Position',[100 100 900 900])
set(gcf, 'PaperPosition', [2 1 15 10]);
print('-depsc2','-r0',filename);

% Read the trajectory file
tra = read_trajectory('trajectory.lsl');

% Open a new figure and set the geographical projection and region
figure(1);
clf;
load coast
h=axesm('MapProjection','stereo','origin',[ 90 70 ]);
gridm;
h=plotm(lat,long,'Color','k','LineWidth',1.5)
%axis([-1.5 1.5 -2 0]);

% Select a set of trajectories
tra.select = (tra.lon > 6) & (tra.lon < 10) & (tra.lat > 44) & (tra.lat < 49) & (tra.time == 0);
tra.select = tra.label(tra.select);


% Plot the trajectories
col.shift       = 47;
col.name        = [ pwd '/col_table.txt' ];
col.spacing     = 900:-50:200;
col.field       = 'p';
col.orientation = 'vert';
plot_trajectory(tra,col,[0 -48 -96]);

% Save figure
filename =   'trajectory.eps';
set(gcf, 'PaperPosition', [2 1 15 10]);
print('-depsc2','-r0',filename);
