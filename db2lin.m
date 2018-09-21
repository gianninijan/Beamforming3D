function [y] = db2lin(x)
 %db conversion to linear scale

    y = 10.^(x./10);
    
end