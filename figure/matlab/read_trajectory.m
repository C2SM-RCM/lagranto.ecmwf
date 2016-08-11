% ----------------------------------------------------------------------
% Read a trajectory file
% ----------------------------------------------------------------------

function out = read_trajectory(filename);

% Open the trajectory file
fid = fopen(filename,'r')

% Read the info line (tra.info)
line     = textscan(fid,'%s',1,'delimiter','\n');
tra.info = line{1};

% Read the line with field description (skip empty lines)
textscan(fid,'%s',1,'delimiter','\n');
line     = textscan(fid,'%s',1,'delimiter','\n');
textscan(fid,'%s',1,'delimiter','\n');
line     = char(line{1});

% Get the names of the fields (tra.field, tra.nfield)
indl       = 1;
indr       = 1;
tra.nfield = 0;
while ( indl <= length(line) )
    while ( line(indl) == ' '  )
       indl=indl+1;
    end
    indr=indl;
    while ( (line(indr) ~= ' ') & (indr < length(line)) )
       indr=indr+1;
    end
    tra.nfield            = tra.nfield + 1;
    tra.field(tra.nfield) = cellstr(line(indl:indr));
    indl = indr + 1;
end

% Get the number of times for a single trajectory (tra.ntime)
fpos = ftell(fid);
textscan(fid,'%s',1,'delimiter','\n');
line = textscan(fid,'%s',1,'delimiter','\n');
tra.ntime = 0;
while ( ~ strcmp(cellstr(line{1}),'') )
    line = textscan(fid,'%s',1,'delimiter','\n');   
    tra.ntime = tra.ntime + 1;
end    
fseek(fid,fpos,'bof');

% Read the trajectories
format = '';
for i=1:tra.nfield
   format = [ format '%n' ];
end
lines = textscan(fid,format);

% Set the total number of trajectories
tra.ntra = length(lines{1})/tra.ntime;

% Attribute the field names to the columns
for i=1:tra.nfield
    name = char(tra.field(i));
    tra.(name) = lines{i};
end

% Attribute a unique label to each trajectory
tra.label = 1+floor( (0:tra.ntra*tra.ntime-1) / tra.ntime )'; 

% Close trajectory file
fclose(fid);

% Return the output structure
out = tra;


