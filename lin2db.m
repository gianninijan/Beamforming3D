function [y] = lin2db(x)
 %conversion to linear scale dB

    y = 10*log10(x);
    
end