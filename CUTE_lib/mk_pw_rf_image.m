function[I] = mk_pw_rf_image(xvec , zvec , c , phi , rf_signal , ...
    t_array , axiscor , envelope_on)
if nargin == 7
    envelope_on = 1;
end
[xs , zs] = meshgrid(xvec , zvec);
Imsize = size(xs);
xs = xs(:);
zs = zs(:);
[N_elements , ~] = size(rf_signal);
delayts = (zs*cos(phi)+xs*sin(phi)) / c;
I = zeros(length(xs) , 1);
for elecount = 1 : N_elements
    linetemp = rf_signal(elecount , :);
    if(envelope_on)
        linetemp = hilbert(linetemp);
    end
    delayrs = sqrt((xs - axiscor(elecount)).^2 + zs.^2) / c;
    interptemp = interp1(t_array , linetemp , delayts + ...
    delayrs);
    interptemp(isnan(interptemp)) = 0;
    apo_fnum = (zs ./ abs(axiscor(elecount) - xs)) > 0.5;
    interptemp = interptemp .* (1-isnan(interptemp));
    I = I + interptemp .* apo_fnum;
%     I = I + interptemp;
end
I = reshape(I , Imsize);
end