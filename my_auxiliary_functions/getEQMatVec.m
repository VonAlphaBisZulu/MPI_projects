% 
% find(iJO1366.stoichMat(strcmp(cellstr(iJO1366.specID),'succ_c'),:))
function [ output_args ] = getEQ2( A, m, q, reacnr, spec )
if exist('spec','var')
    if spec
        if ischar(reacnr)
            reacnr = find(A(strcmp(m,reacnr),:));
        elseif ismatrix(reacnr)
            reacnr = find(A(reacnr,:));
        end
    end
end
if ischar(reacnr)
    reacnr = find(strcmp(q,reacnr),1);
end
for i = 1:length(reacnr)
    zw=find(A(:,reacnr(i))<0);
    if(~isempty(zw))
            str=[num2str(-A(zw(1),reacnr(i))),' ',deblank(m(zw(1),:))];
            for j=2:length(zw)
                    str=[str,' + ',num2str(-A(zw(j),reacnr(i))),' ',deblank(m(zw(j),:))];
            end
    else
        str = '';
    end            

    if(~strcmp('mue',deblank(q(reacnr(i)))))
        str=[str,' = '];
    end
    zw=find(A(:,reacnr(i))>0);
    if(~isempty(zw))
            str=[str,num2str(A(zw(1),reacnr(i))),' ',deblank(m(zw(1),:))];
            for j=2:length(zw)
                str=[str,' + ',num2str(A(zw(j),reacnr(i))),' ',deblank(m(zw(j),:))];
            end                                                  
    end            
    disp(strjoin([strtrim(q(reacnr(i))) '   ' char(9613) '   ' str ]));
end
end

