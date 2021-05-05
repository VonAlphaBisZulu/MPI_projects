function numcode = getSlurmArrayCode(product,coupling,maxCost,gene_mcs,atpm)
% digits 1-4  : product
% digits 5-7 : weak coupling, directional coupling, substrate coupling
% digits 8-16: maxCost
% digits 17: gene MCS
% digits 18: ATPM
    switch coupling
        case 'auto'
            numcode_coupling = '000';
        case 'potential'
            numcode_coupling = '001';
        case 'weak'
            numcode_coupling = '010';
        case 'directional'
            numcode_coupling = '011';
        case 'ATP'
            numcode_coupling = '100';
        case 'substrate'
            numcode_coupling = '101';
    end
    numcode_product = dec2bin(product,4);
    numcode_maxCost = dec2bin(maxCost,9);
    numcode = [numcode_product numcode_coupling numcode_maxCost num2str(gene_mcs) num2str(atpm)];
    numcode = bin2dec(numcode);
end

