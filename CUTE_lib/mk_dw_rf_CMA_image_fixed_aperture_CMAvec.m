function[I , mask] = mk_dw_rf_CMA_image_fixed_aperture_CMAvec(imxvec , imzvec , c , emitx ,...
    rf_signal , t_array , axiscor , CMAvec , envelope_on)
if nargin == 8
    envelope_on = 1;
end
apersize = axiscor(end) - axiscor(1);

[N_elements , ~] = size(rf_signal);
[imxs , imzs] = meshgrid(imxvec , imzvec); Imsize = size(imxs);
[nz , nx] = size(imxs);

tof = zeros(N_elements , nz , nx);
[~ , emitxn] = min(abs(axiscor - emitx));
for elecount = 1 : N_elements
    tof(elecount , : , :) = sqrt((imxs - axiscor(elecount)).^2 + imzs.^2) ./ c;
end
cc=1;
for CMA = CMAvec
    emitangles = atan((imxs-emitx)./imzs);
    recvangles = 2 * CMA - emitangles;
    recvangles(recvangles <= - pi/2) = - 49/100 * pi;
    recvangles(recvangles >=  pi/2) =  49/100 * pi;
    
    maxapersize = 0.4 * apersize;
    xmin = imxs - imzs .* tan(recvangles) - maxapersize / 2;
    xmax = imxs - imzs .* tan(recvangles) + maxapersize / 2;
    xmin = max(axiscor(1) , xmin); xmax = min(axiscor(end) , xmax);

    Itemp = imxs * 0;
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
        Itemp = Itemp + interptemp .* apo .* recvhann .* apo_fnum;
        totalapo = totalapo | apo;
    end
    masktemp = reshape(totalapo , Imsize);
    I(cc , : , :) = Itemp;
    mask(cc , : , :) = masktemp;
    cc = cc + 1;
end
end