function TF = contains(str,pattern)
if ~iscell(pattern)
    pattern = cellstr(pattern);
end
if ~iscell(str)
    str = cellstr(str);
end
if all(size(pattern)==[1 1])
    TF = ~cellfun(@isempty,strfind(str,char(pattern)));
else
    TF = false(size(str));
    for idx = 1:numel(pattern)
        TF = or(TF , strcmp(str,pattern(idx)));
    end
end