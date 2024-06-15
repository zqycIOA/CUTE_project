function[I , mask] = mk_dw_rf_CMA_image_fixed_aperture(xvec , zvec , c , emitx , rf_signal , ...
    t_array , axiscor , CMA , envelope_on)
if nargin == 8
    envelope_on = 1;
end
apersize = axiscor(end) - axiscor(1);

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
maxapersize = 0.1 * apersize;
xmin = xs - zs .* tan(recvangles) - maxapersize / 2;
xmax = xs - zs .* tan(recvangles) + maxapersize / 2;
xmin = max(axiscor(1) , xmin); xmax = min(axiscor(end) , xmax);
for elecount = 1 : N_elements
    delayrs(:,elecount)=sqrt((xs - axiscor(elecount)).^2 + zs.^2) / c;
end
I = zeros(length(xs) , 1);
totalapo = xs * 0;
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
%     I = I + interptemp .* apo .* recapo(elecount);
    I = I + interptemp .* apo.*recvhann.*apo_fnum;
    totalapo = totalapo | apo;
end
% I = hilbert(reshape(I , Imsize));
I = reshape(I , Imsize);
mask = reshape(totalapo , Imsize);
end