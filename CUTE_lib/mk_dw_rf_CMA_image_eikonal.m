function[I , mask] = mk_dw_rf_CMA_image_eikonal(imxvec , imzvec , cxvec , czvec , Cmap , emitx ,...
    rf_signal , t_array , axiscor , CMA , envelope_on)
if nargin == 10
    envelope_on = 1;
end
bfgap = 0.1e-3;
[N_elements , ~] = size(rf_signal);
[imxs , imzs] = meshgrid(imxvec , imzvec); Imsize = size(imxs);
[cxs , czs] = meshgrid(cxvec , czvec);
imxvec_bf = imxvec(1) : bfgap : imxvec(end);
imzvec_bf = 0 : bfgap : imzvec(end);
[xs_bf , zs_bf] = meshgrid(imxvec_bf , imzvec_bf);
[~ , emitxn] = min(abs(axiscor - emitx));
Cmap_bf = interp2_padded(cxs , czs , Cmap , xs_bf , zs_bf);
% Cmap_bf = interp2(cxs , czs , Cmap , xs_bf , zs_bf);
% Cmap_bf(isnan(Cmap_bf)) = median(median(Cmap(1,:)));
for elecount = 1 : N_elements
    elex = axiscor(elecount);
    [~ , devn] = min(abs(elex - imxvec_bf));
    cxdev = elex - imxvec_bf(devn);
    Cmap_dev = interp2(xs_bf, zs_bf, Cmap_bf, xs_bf + cxdev, zs_bf, 'spline');
    toftemp = bfgap * msfm2d(Cmap_dev, [1; devn], true, true); 
    tof(elecount , : , :) = interp2(xs_bf + cxdev, zs_bf, toftemp, imxs, imzs, 'spline');
    tof1(elecount , : , :) = sqrt((imxs - axiscor(elecount)).^2 + imzs.^2) ./ Cmap(1);
end

emitangles = atan((imxs-emitx)./imzs);
recvangles = 2 * CMA - emitangles;
recvangles(recvangles <= - pi/2) = - 49/100 * pi;
recvangles(recvangles >=  pi/2) =  49/100 * pi;
pitch = axiscor(2) - axiscor(1); maxapersize = 50 * pitch;
xmin = imxs - imzs .* tan(recvangles) - maxapersize / 2;
xmax = imxs - imzs .* tan(recvangles) + maxapersize / 2;
xmin = max(axiscor(1) , xmin); xmax = min(axiscor(end) , xmax);

I = imxs * 0;
totalapo = imxs * 0;
for elecount = 1 : N_elements
    toftemp = squeeze(tof(elecount,:,:) + tof(emitxn,:,:));
    linetemp = rf_signal(elecount , :);
    if(envelope_on)
        linetemp = hilbert(linetemp);
    end
    realapersize = xmax-xmin;
    recvhann = 0.5 + 0.5.*cos(((xmin + xmax)./2 - axiscor(elecount))./...
        maxapersize * pi * 2);
    apo = (axiscor(elecount) > xmin) & (axiscor(elecount) < xmax)&...
        (realapersize/maxapersize > 0.8);
    apo_fnum = (imzs ./ abs(axiscor(elecount) - imxs)) > 1;
    
    interptemp = interp1(t_array , linetemp , toftemp);
    interptemp = interptemp .* (1-isnan(interptemp));
    I = I + interptemp .* apo .* recvhann .* apo_fnum;
    totalapo = totalapo | apo;
end
mask = reshape(totalapo , Imsize);
end