function[I , mask] = mk_pw_rf_CMA_image_CMAvec(xvec , zvec , c , phi , rf_signal , ...
    t_array , axiscor , CMAvec , envelope_on)
if nargin == 8
    envelope_on = 1;
end
[xs , zs] = meshgrid(xvec , zvec);
Imsize = size(xs);
xs = xs(:);
zs = zs(:);
[N_elements , ~] = size(rf_signal);
delayts = (zs*cos(phi)+xs*sin(phi)) / c;
I0 = zeros(length(xs) , 1);
for elecount = 1 : N_elements
    linetemp = rf_signal(elecount , :);
    if(envelope_on)
        linetemp = hilbert(linetemp);
    end
    delayrs = sqrt((xs - axiscor(elecount)).^2 + zs.^2) / c;
    interptemp = interp1(t_array , linetemp , delayts + ...
    delayrs);
    interptemp(isnan(interptemp)) = 0;
    apo_fnum = (zs ./ abs(axiscor(elecount) - xs)) > 0.8;
    interptemp = interptemp .* (1-isnan(interptemp));
    I0 = I0 + interptemp .* apo_fnum;
%     I0 = I0 + interptemp;
end
I0 = reshape(I0 , Imsize);
cc = 1;
for CMA = CMAvec
     Itemp = k_space_filt(I0 , xvec , zvec , CMA);
%      Itemp = I0;
     I(cc , : , :) = Itemp;
     masktemp = (xs <= min(axiscor(end) , axiscor(end) - zs.* tan(CMA))) &...
         (xs >= max(axiscor(1) , axiscor(1) - zs.* tan(CMA)));
     mask(cc , : , :) = reshape(masktemp , Imsize);
     cc = cc + 1;
end
end