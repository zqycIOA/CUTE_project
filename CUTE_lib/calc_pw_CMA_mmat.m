function[Mmat] = calc_pw_CMA_mmat(totalmask , ROIparam , psparam , arrayparam)
    xvec = ROIparam.xvec; zvec = ROIparam.zvec;
    [xs , zs] = meshgrid(xvec , zvec);
    sigmalength = length(xs(:));
    CMAvec = psparam.CMAvec;
    
    sigmavec = 1 : sigmalength;
    matlength = 0;
    phivec = psparam.psstart:psparam.psgap:psparam.psend;
    emitsize = length(phivec);
    for pc = 1 : length(totalmask(:))
        matlength = matlength + sum(sum(totalmask{pc}));
    end
    Mmat = zeros(matlength, sigmalength);
    cc = 1; mc = 1;
    for CMA = CMAvec
        for emitcount = 1:emitsize-1
            masktemp = totalmask{cc , emitcount};
            for pixcount = sigmavec(masktemp == 1)
                recvangle1 = 2 * CMA - phivec(emitcount);
                recvangle2 = 2 * CMA - phivec(emitcount + 1);
                wvec1 = calc_pw_w_vec(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    phivec(emitcount)) + ...
                    calc_pw_w_vec(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    recvangle1);
                wvec2 = calc_pw_w_vec(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    phivec(emitcount + 1)) + ...
                    calc_pw_w_vec(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    recvangle2);
                Mmat(mc , :) = wvec1 - wvec2;
                mc = mc + 1;
            end
        end
        cc = cc + 1;
    end
end