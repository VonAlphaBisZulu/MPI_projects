load('C:\Users\Philipp\Documents\Python\mcs_python.mat')
name_mcs = cellstr(whatever_data);
mcs_python = zeros(gcnap.numr,numel(name_mcs));
for i = 1:numel(name_mcs)
    kos = strsplit(strtrim(whatever_data(i,:)),',');
    for j = 1:numel(kos)
        if ~strcmp(kos{j},'EX_o2_e')
            kos{j} = ['GP-' kos{j}];
        end
    end
    mcs_python(ismember(cellstr(gcnap.reacID),kos),i) = -1;
end

%% Verify MCS
% target module
modules{1}.type = 'lin_constraints';
modules{1}.sense = 'desired';
modules{1}.V = gcnap.mcs.D{:};
modules{1}.v = gcnap.mcs.d{:};
modules{2}.type = 'lin_constraints';
modules{2}.sense = 'target';
modules{2}.V = gcnap.mcs.T{:};
modules{2}.v = gcnap.mcs.t{:};

valid = verify_mcs(gcnap,modules,mcs_python,3,-1);