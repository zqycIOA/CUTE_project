function [totalps] = calc_pw_ps_prediction(ROIparam , psparam , arrayparam) 
    xvec = ROIparam.xvec; zvec = ROIparam.zvec;
    [xs , zs] = meshgrid(xvec , zvec);
    Imsize = size(xs);
    sigmalength = length(xs(:));
    Cmap0 = ROIparam.Cmap0; Cmap0 = Cmap0(:);
    c = ROIparam.c;
    ds = 1 ./ c - 1 ./ Cmap0;
    axiscor = arrayparam.axiscortx;
    phivec = psparam.psstart:psparam.psgap:psparam.psend;
    emitsize = length(phivec);
    if isfield(psparam , 'CMAvec')
        CMAvec = psparam.CMAvec;
        cc = 1;
        for CMA = CMAvec
            ec = 1;
            for emitcount = 1:emitsize-1
                Mmattemp = zeros(sigmalength, sigmalength);
                mc = 1;
                recvangle1 = 2 * CMA - phivec(emitcount);
                recvangle2 = 2 * CMA - phivec(emitcount + 1);
                for pixcount = 1 : sigmalength
                    wvec1 = calc_pw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                        phivec(emitcount) , 0) + ...
                        calc_pw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                        recvangle1 , 0);
                    wvec1 = wvec1(:);
                    wvec2 = calc_pw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                        phivec(emitcount + 1) , 0) + ...
                        calc_pw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                        recvangle2 , 0);
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
                wvec1 = calc_pw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    phivec(emitcount) , 0);
                wvec1 = wvec1(:);
                wvec2 = calc_pw_w_mat(xvec , zvec , xs(pixcount) , zs(pixcount),...
                    phivec(emitcount + 1) , 0);
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