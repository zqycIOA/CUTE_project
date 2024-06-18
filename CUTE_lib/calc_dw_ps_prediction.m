function [totalps] = calc_dw_ps_prediction(ROIparam , psparam , arrayparam) 
    xvec = ROIparam.xvec; zvec = ROIparam.zvec;
    [xs , zs] = meshgrid(xvec , zvec);
    Imsize = size(xs);
    sigmalength = length(xs(:));
    Cmap0 = ROIparam.Cmap0; Cmap0 = Cmap0(:);
    c = ROIparam.c;
    ds = 1 ./ c - 1 ./ Cmap0;
    axiscor = arrayparam.axiscortx;
    emitvec = axiscor(psparam.psstart:psparam.psgap:psparam.psend);
    emitsize = length(emitvec);
    if isfield(psparam , 'CMAvec')
        CMAvec = psparam.CMAvec;
        cc = 1;
        for CMA = CMAvec
            ec = 1;
            for emitcount = 1:emitsize-1
                Mmattemp = zeros(sigmalength, sigmalength);
                mc = 1;
                for pixcount = 1 : sigmalength
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
                    Mmattemp(mc , :) = wvec1 - wvec2;
                    mc = mc + 1;
                end
                pstemp = Mmattemp * ds;
                pstemp = reshape(pstemp , Imsize);
                totalps{cc , ec} = pstemp;
                ec = ec + 1;
            end
            cc = cc + 1;
        end
    else
        ec = 1;
        for emitcount = 1:emitsize-1
            Mmattemp = zeros(sigmalength, sigmalength);
            mc = 1;
            for pixcount = 1 : sigmalength
                wvec1 = calc_dw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    emitvec(emitcount) , 0);
                wvec1 = wvec1(:);
                wvec2 = calc_dw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    emitvec(emitcount + 1) , 0);
                wvec2 = wvec2(:);
                Mmattemp(mc , :) = wvec1 - wvec2;
                mc = mc + 1;
            end
            pstemp = Mmattemp * ds;
            pstemp = reshape(pstemp , Imsize);
            totalps{ec} = pstemp;    
            ec = ec + 1;
        end
    end
end