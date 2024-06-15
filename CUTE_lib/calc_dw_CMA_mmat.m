function[Mmat] = calc_dw_CMA_mmat(totalmask , ROIparam , psparam , arrayparam)
    xvec = ROIparam.xvec; zvec = ROIparam.zvec;
    [xs , zs] = meshgrid(xvec , zvec);
    sigmalength = length(xs(:));
    CMAvec = psparam.CMAvec;
    
    sigmavec = 1 : sigmalength;
    matlength = 0;
    axiscor = arrayparam.axiscortx;
    emitvec = axiscor(psparam.psstart:psparam.psgap:psparam.psend);
    emitsize = length(emitvec);
    for pc = 1 : length(totalmask(:))
        matlength = matlength + sum(sum(totalmask{pc}));
    end
    Mmat = zeros(matlength, sigmalength);
    cc = 1; mc = 1;
    for CMA = CMAvec
        for emitcount = 1:emitsize-1
            masktemp = totalmask{cc , emitcount};
            for pixcount = sigmavec(masktemp == 1)
                emitangle1 = atan((xs(pixcount) - emitvec(emitcount))/zs(pixcount));
                emitangle2 = atan((xs(pixcount) - emitvec(emitcount + 1))/zs(pixcount));
                recvangle1 = 2 * CMA - emitangle1;
                recvangle2 = 2 * CMA - emitangle2;
                recvangle1(recvangle1 <= - pi/2) = - 49/100 * pi;
                recvangle2(recvangle2 >=  pi/2) =  49/100 * pi;
                recvpos1 = xs(pixcount) - zs(pixcount) * tan(recvangle1);
                recvpos2 = xs(pixcount) - zs(pixcount) * tan(recvangle2);
                wvec1 = calc_dw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    emitvec(emitcount) , 0) + ...
                    calc_dw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    recvpos1 , 0);
                wvec1 = wvec1(:);
                wvec2 = calc_dw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    emitvec(emitcount + 1) , 0) + ...
                    calc_dw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    recvpos2 , 0);
                wvec2 = wvec2(:);
                Mmat(mc , :) = wvec1 - wvec2;
                mc = mc + 1;
            end
        end
        cc = cc + 1;
    end
end