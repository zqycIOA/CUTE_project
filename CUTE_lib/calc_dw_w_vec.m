function[w] = calc_dw_w_vec(xvec, zvec, x, z, elex)

nx = length(xvec);
nz = length(zvec);
dx = xvec(2)-xvec(1);
dz = zvec(2)-zvec(1);

xstart = xvec(1);
zstart = zvec(1);
xend = xvec(end);
zend = zvec(end);
%% Define & Plot Target Position and Ray

i2 = round((x - xstart) / dx) + 1;
j2 = round((z - zstart) / dz) + 1;
phi = acot(z/(x-elex));
x1 = elex + zstart*tan(phi);
z1 = zstart;
%% Calculate Intersection with x and z Planes
% x Plane
alphaxs = ((1:j2-1)-1)/(j2-1);
xxps = x + (alphaxs - 1) .* tan(phi) * (z - z1);
zxps = alphaxs .* z + (1 - alphaxs) .* zstart;
% scatter(xxps , zxps ,'r','x' );
xflag = alphaxs * 0 + 1;

% z Plane
Deltax1 = abs(tan(phi) * abs(z - z1)) - floor(abs(tan(phi) * (z - z1)) / dx) * dx;
kmax = floor(abs((z - z1) * tan(phi) / dx));
alphazs = (((1 : kmax) - 1) * dx + Deltax1) / abs(tan(phi) * (z - z1));
xzps = x + (alphazs - 1) .* tan(phi) * (z - z1);
zzps = alphazs .* z + (1 - alphazs) .* zstart;
% scatter(xzps , zzps ,'g' ,'x' );
zflag = alphazs * 0 + 2;

%% Combine Intersections and Trim

totalxs = [xxps xzps];
totalzs = [zxps zzps];
totalalphas = [alphaxs alphazs];
totalflag = [xflag zflag];
[totalalphas , I] = sort(totalalphas);
totalxs = totalxs(I);
totalzs = totalzs(I);
totalflag = totalflag(I);

trimind = (totalxs >= xstart) & (totalxs <= xend);
totalxs = totalxs(trimind);
totalzs = totalzs(trimind);
totalalphas = totalalphas(trimind);
totalflag = totalflag(trimind);

if ~isempty(totalalphas)
    totalalphas = (totalalphas - totalalphas(1)) ./ (1 - totalalphas(1));

    crossnum = length(totalxs);

    totalsteps = ([totalalphas(2 : end) 1] - totalalphas) * sqrt(1 + tan(phi)^2)...
        * (z - totalzs(1));
    totalsteps(1) = totalsteps(1) + sqrt((totalxs(1) - elex).^2 + totalzs(1).^2);
end


%% Calculate w Matrix
w = zeros(nz , nx);
if ~isempty(totalalphas)
    edgeflag = ((abs(totalxs - xstart) < 1e-5) | ((abs(totalxs - xend)) < 1e-5));
    for crosscount = 1 : crossnum
        if (totalflag(crosscount) == 1)
            if ~(edgeflag(crosscount))
                xind_low = floor((totalxs(crosscount) - xstart) / dx) + 1;
                xpos_low = xstart + dx * (xind_low - 1);
                xind_high = xind_low + 1;
                xpos_high = xstart + dx * (xind_high - 1);
                zind = round((totalzs(crosscount) - zstart) / dz) + 1;
                w(zind , xind_high) = w(zind , xind_high) + (totalxs(crosscount)...
                    - xpos_low)/dx * totalsteps(crosscount);
                w(zind , xind_low) = w(zind , xind_low) + (xpos_high - ...
                    totalxs(crosscount))/dx * totalsteps(crosscount);
            else
                 xind = round((totalxs(crosscount) - xstart) / dx) + 1;
                 zind = round((totalzs(crosscount) - zstart) / dz) + 1;
                 w(zind , xind) = totalsteps(crosscount);
            end
        elseif (totalflag(crosscount) == 2)
            zind_low = floor((totalzs(crosscount) - zstart) / dz) + 1;
            zpos_low = zstart + dz * (zind_low - 1);
            zind_high = zind_low + 1;
            zpos_high = zstart + dz * (zind_high - 1);
            xind = round((totalxs(crosscount) - xstart) / dx) + 1;
            w(zind_high , xind) = w(zind_high , xind) + (totalzs(crosscount)...
                - zpos_low)/dz * totalsteps(crosscount);
            w(zind_low , xind) = w(zind_low , xind) + (zpos_high - ...
                totalzs(crosscount))/dz * totalsteps(crosscount);
        end
    end
else
    w(j2 , i2) = sqrt((xvec(i2) - elex).^2 + zvec(j2).^2);
end

w = w(:);


end