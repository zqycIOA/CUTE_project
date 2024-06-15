function[V] = interp2_padded(vx , vz , v , x , z)
    V_linear = interp2(vx , vz , v , x , z,'linear');
    V_nanind = isnan(V_linear);
    V_linear(V_nanind) = 0;
    V1st = mean(v(1:5,:));
    vbound = interp1(vx(1 , :) , V1st , x(1 , :) , 'spline');
%     vbound = vbound * 0 + mean(vbound);
    V_nearest = (x(:,1) * 0 + 1) * vbound;
    V = V_linear + V_nearest .* V_nanind;
end