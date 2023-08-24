function [xc,yc,zc] = candidate_voxel(v,w,data_size)
    x = v(1);
    y = v(2);
    z = v(3);
    
    xt = [x-floor(w(1)/2):x-1 x x+1:x+floor(w(1)/2)];
    yt = [y-floor(w(2)/2):y-1 y y+1:y+floor(w(2)/2)];
    zt = [z-floor(w(3)/2):z-1 z z+1:z+floor(w(3)/2)];
    
    if x == 1
        xt = [x x+1:x+floor(w(1)/2)];
    elseif x-floor(w(1)/2) <= 0 && x > 1
        xt = [1:x-1 x x+1:x+floor(w(1)/2)];
    elseif x+floor(w(1)/2) > data_size(1)
        xt = [x-floor(w(1)/2):x-1 x x+1:data_size(1)];
    elseif x == data_size(1)
        xt = [x-floor(w(1)/2):x-1 x];
    end
    
    if y == 1
        yt = [y y+1:y+floor(w(2)/2)];
    elseif y-floor(w(2)/2) <= 0 && y > 1
        yt = [1:y-1 y y+1:y+floor(w(2)/2)];
    elseif y+floor(w(2)/2) > data_size(2)
        yt = [y-floor(w(2)/2):y-1 y y+1:data_size(2)];
    elseif y == data_size(2)
        yt = [y-floor(w(2)/2):y-1 y];
    end
    
    if z == 1
        zt = [z z+1:z+floor(w(3)/2)];
    elseif z-floor(w(3)/2) <= 0 && z > 1
        zt = [1:z-1 z z+1:z+floor(w(3)/2)];
    elseif z+floor(w(3)/2) > data_size(3)
        zt = [z-floor(w(3)/2):z-1 z z+1:data_size(3)];
    elseif z == data_size(3)
        zt = [z-floor(w(3)/2):z-1 z];
    end
    
    
    [xt2,yt2,zt2] = meshgrid(xt,yt,zt);
    
    xc = xt2(:);
    yc = yt2(:);
    zc = zt2(:);
    
    [m,~] = size(xc);
    
    for i = 1:m
        if (xc(i)==v(1)) && (yc(i)==v(2)) && (zc(i)==v(3))
            xc(i) = [];
            yc(i) = [];
            zc(i) = [];
            
            break;
        end
    end
    
    