function [y] = lin2dbm(x)
 %conversion to linear scale dBm

    y = 10*log10(x./1e-3);
    
end