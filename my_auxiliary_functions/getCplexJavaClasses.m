cpo= javaObject('ilog.cplex.IloCplex');
ilo_cplex_class= cpo.getClass();
%A# java.lang.Class.forName('ilog.cplex.IloCplex') throws ClassNotFoundException ?!?
cl= ilo_cplex_class.getClasses();
sn= cell(size(cl));
for i= 1:length(cl)
sn{i}= char(cl(i).getSimpleName());
end
    
while isstruct(getfield(ci, paramTable{i(1),1:i(2)})) && (i(1) <= size(sn,1))
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