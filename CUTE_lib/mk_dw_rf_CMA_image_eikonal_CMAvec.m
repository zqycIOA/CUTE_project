function[I , mask] = mk_dw_rf_CMA_image_eikonal_CMAvec(xvec , zvec , c , emitx , rf_signal , ...
    t_array , axiscor , CMA , envelope_on)
if nargin == 8
    envelope_on = 1;
end
pitch = axiscor(2) - axiscor(1);

[xs , zs] = meshgrid(xvec , zvec);
Imsize = size(xs);
xs = xs(:);
zs = zs(:);
[N_elements , ~] = size(rf_signal);
delayts = sqrt((xs - emitx).^2 + zs.^2) / c;
delayrs = zeros(length(zs) , N_elements);
emitangles = atan((xs-emitx)./zs);
recvangles = 2 * CMA - emitangles;
recvangles(recvangles <= - pi/2) = - 49/100 * pi;
recvangles(recvangles >=  pi/2) =  49/100 * pi;
maxapersize = 100 * pitch;
xmin = xs - zs .* tan(recvangles) - maxapersize / 2;
xmax = xs - zs .* tan(recvangles) + maxapersize / 2;
xmin = max(axiscor(1) , xmin); xmax = min(axiscor(end) , xmax);
for elecount = 1 : N_elements
    delayrs(:,elecount)=sqrt((xs - axiscor(elecount)).^2 + zs.^2) / c;
end

cc = 1;
    for CMA = CMAvec
        emitangles = atan((imxs-emitx)./imzs);
        recvangles = 2 * CMA - emitangles;
        recvangles(recvangles <= - pi/2) = - 49/100 * pi;
        recvangles(recvangles >=  pi/2) =  49/100 * pi;
        pitch = axiscor(2) - axiscor(1); maxapersize = 100 * pitch;
        xmin = imxs - imzs .* tan(recvangles) - maxapersize / 2;
        xmax = imxs - imzs .* tan(recvangles) + maxapersize / 2;
        xmin = max(axiscor(1) , xmin); xmax = min(axiscor(end) , xmax);
    
        Itemp = imxs * 0;
        totalapo = imxs * 0;
        for elecount = 1 : N_elements
        %     apo = (axiscor(elecount) > xmin) & (axiscor(elecount) < xmax);
            realapersize = xmax-xmin;
            recvhann = 0.5 + 0.5.*cos(((xmin + xmax)./2 - axiscor(elecount))./...
                maxapersize * pi * 2);
            apo = (axiscor(elecount) > xmin) & (axiscor(elecount) < xmax)&...
                (realapersize/maxapersize > 0.8);
        %     apo = 1;
            linetemp = rf_signal(elecount , :);
            if(envelope_on)
                linetemp = hilbert(linetemp);
            end
            interptemp = interp1(t_array , linetemp , delayts + ...
            delayrs(: , elecount));
            apo_fnum = (zs ./ abs(axiscor(elecount) - xs)) > 1;
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
