function f = roundDec( f , decimals )
    switch nargin
        case 1
            f = round(f);
        case 2
            f = round(f*(10^(decimals)))/(10^(decimals));
        otherwise
            error('Invalid number of arguments');
    end
end

