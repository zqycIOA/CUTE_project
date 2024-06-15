function[I] = mk_dw_rf_image(xvec , zvec , c , emitx , rf_signal , ...
    t_array , axiscor , envelope_on)
if nargin == 7
    envelope_on = 1;
end
fnum = 0.5;
% recapo = kaiser(length(axiscor) , 2) * 0  + 1 ;
[xs , zs] = meshgrid(xvec , zvec);
Imsize = size(xs);
xs = xs(:);
zs = zs(:);
[N_elements , ~] = size(rf_signal);
delayts = sqrt((xs - emitx).^2 + zs.^2) / c;
I = zeros(length(xs) , 1);
% apo = tukeywin(N_elements,0.4);
for elecount = 1 : N_elements
%     apo = 1;
    linetemp = rf_signal(elecount , :);
    if(envelope_on)
        linetemp = hilbert(linetemp);
    end
    delayrs = sqrt((xs - axiscor(elecount)).^2 + zs.^2) / c;
    interptemp = interp1(t_array , linetemp , delayts + ...
    delayrs);
%     apo_fnum = (zs ./ abs(axiscor(elecount) - xs)) > fnum;
    interptemp(isnan(interptemp)) = 0;
%     I = I + interptemp .* apo(elecount);
    I = I + interptemp;
end
I = reshape(I , Imsize);
end