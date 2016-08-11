function output = col_load(coltable, spacing, shift, varargin)
% Routine handles colortable
%
% coltable  input    string contains ctb-file and its path
%                    e.g. /home/user/matlab/.../file
%
%                    if only file: path = /home/mischa/matlab/ctb/
% 
% date2     input    shifter for colortable
% 
% varargin  input    pointer to define a white space                 
%                  
% output    output   nothing
% 
% ----------------------------------------------------------
%                            Mischa Croci Maspoli (Oct 2004)
% ----------------------------------------------------------

%check if varargin = 1
if length(varargin) == 1
 white=varargin{1};
end
if length(varargin) == 0
 white=0;
end


ctbdefault = '/usr/local/matlabtools/ive_ct/';


% read colortable
if coltable(1) == '/'
 ivecol = load(coltable,'-ascii');
else
 ivecol = load([ctbdefault coltable],'-ascii');
end

% devide values by 255
ivecol = 1./255.*ivecol;

%ivecol=flipdim(ivecol,1);

% number of colors
ncol=length(spacing)-2;

% shifter
shcol=shift;

% define white space
if  length(varargin) == 1
  newmap = ivecol(shcol:ncol+shcol,:);
  newmap(white,:) = 1.;
else
  newmap = ivecol(shcol:ncol+shcol,:);
end

%data=[1:ncol+1;1:1:ncol+1]';
colormap(newmap(:,:));
