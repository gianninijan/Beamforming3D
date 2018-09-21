function [y] = dbm2lin(x)
 %dbm conversion to linear scale

    y = 10.^(x./10 - 3);
    
end