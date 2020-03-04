function [ parts ] = split( str,delimiter )
% splits char or cell into multiple parts
% returns them as cell array
% str: original string
% delimiter: character (array) where to split (default is whitespace ' ')

% prepare parameters
if nargin == 1
    delimiter = ' ';
end

if ischar(delimiter)
    
elseif iscellstr(delimiter)
    delimiter = char(delimiter);
else
    error('split: invalid type of delimiter. Enter char or cellstr');
end

if ischar(str)
    
elseif iscellstr(str)
    str = char(str);
else
    error('split: invalid type of input string. Enter char or cellstr');
end
if isempty(str)
    parts = {''};
    return;
end

% replace multiple delimiters placed in line.
while ~isempty(strfind(str, [delimiter delimiter]))
    str = strrep(str,[delimiter delimiter],delimiter);
end
% start splitting
pos = strfind(str,delimiter)';
pos      = [1;pos+length(delimiter)];
pos(:,2) = [pos(2:end)-length(delimiter)-1;length(str)];
for i=1:size(pos,1)
    parts(i) = {str(pos(i,1):pos(i,2))};
end
parts = parts';
end

