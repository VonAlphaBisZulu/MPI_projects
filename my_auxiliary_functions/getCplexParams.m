function paramTable = getCplexParams( cpx )
ci = setup_cplex_inner_class_access();
paramTable = fieldnames(ci);
i = [1;1];
while isstruct(getfield(ci, paramTable{i(1),1:i(2)})) && (i(1) <= size(paramTable,1))
    fn = fieldnames(ci.(char(paramTable(i(1),i(2)))));
    for j = 0:(length(fn)-1)
        if j>0
            paramTable((i(1)+j):(end+1),i(2)) = paramTable((i(1)+j-1):end,i(2));
        end
        paramTable(i(1)+j,i(2)+1) = fn(j+1);
    end
    i(1) = i(1)+length(fn);
    if i(1) > size(paramTable,1)
        break
    end
end
for j = 1:size(paramTable,1)
    try
        value = cpx.getParam(getfield(ci,paramTable{j,1:2}));
        switch class(value)
            case 'double'
                paramTable(j,3) = num2cell(value);
            case 'logical'
                paramTable(j,3) = num2cell(value);
            case 'char'
                paramTable(j,3) = cellstr(value);
            otherwise
                paramTable(j,3) = cellstr(class(value));
        end
    catch
        disp('not possible')
    end
end