function [model,productID,coupling,maxCost,gene,atpm] = getInfosFormArrayCode(numcode)
    settings = dec2bin(numcode,21);
    switch settings(1)
        case '0'
            model = 'iMLcore';
        case '1'
            model = 'iML1515';
    end
    productID = bin2dec(settings(2:8));
    switch settings([9 10])
        case '00'
            coupling = 'weak growth-coupling';
        case '01'
            coupling = 'directional growth-coupling';
        case '10'
            coupling = 'substrate uptake coupling';
        case '11'
            coupling = 'auto';
    end
    maxCost = bin2dec(settings(11:19));
    gene = bin2dec(settings(20));
    atpm = bin2dec(settings(21));
end
