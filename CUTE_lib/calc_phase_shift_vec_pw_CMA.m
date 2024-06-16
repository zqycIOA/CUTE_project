%------------------------------------------------------------------------%
%psparam: psstart psend psgap pslength CMAvec
% 起始阵元 结束阵元 相差结果间隔阵元 相差提取间隔阵元
%arrayparam: axiscortx axiscorrx
%ROIparam: xvec zvec c
%arrayparam: axiscortx axiscorrx
function varargout = calc_phase_shift_vec_pw_CMA(totalframes , psparam , ROIparam , arrayparam)
    axiscortx = arrayparam.axiscortx;
    axiscorrx = arrayparam.axiscorrx;
    N_lines = length(totalframes);
    phivec = zeros(1 , N_lines);
    for lc = 1 : N_lines
        phivec(lc) = totalframes(lc).phi;
    end
    psstart = psparam.psstart;
    psend = psparam.psend;
    pslength = psparam.pslength;
    xvec = ROIparam.xvec; zvec = ROIparam.zvec; c = ROIparam.c;
    xlength = length(xvec); zlength = length(zvec);

    emitvec = psstart:pslength:psend;
    [~ , emitnvec] = find(abs(phivec - emitvec')<1e-10);
    emitnvec = emitnvec';
    if(isfield(psparam , 'CMAvec'))
    CMAvec = psparam.CMAvec;
    tc = 1;
    for emitn = emitnvec(1 : end - 1)
        [pstemp{tc} , masktemp{tc}] = calc_phase_shift_Hermitian_pw_CMA_CMAvec(xvec ,...
            zvec , c , totalframes(emitnvec(tc)) , totalframes(emitnvec(tc+1)) , axiscorrx , CMAvec);
        tc = tc + 1;
    end
    psgap = psparam.psgap; emitlength = psgap/pslength;
    pstxelen = emitvec(1) : psgap : emitvec(end) - psgap;
    cc = 1;
    for CMA = CMAvec
        ec = 1;
        for txcount = 1:emitlength:length(emitnvec)-1
            ps = zeros(zlength , xlength);
            mask = zeros(zlength , xlength);
            for framecount = txcount : txcount + emitlength-1
                pc = pstemp{framecount};
                mc = masktemp{framecount};
                ps = ps + squeeze(pc(cc , : , :)); 
                mask = mask | squeeze(mc(cc , : , :));
            end
            totalps{cc , ec} = ps;
            totalmask{cc , ec} = double(mask);
            ec = ec + 1
        end
        cc = cc + 1;
    end
    varargout{1} = totalps;
    varargout{2} = totalmask;
    else
    tc = 1;
    for emitn = emitnvec(1 : end - 1)
        pstemp{tc} = calc_phase_shift_Hermitian_pw(xvec ,...
            zvec , c , totalframes(emitnvec(tc)) , totalframes(emitnvec(tc+1)) , axiscorrx);
        tc = tc + 1;
    end
    psgap = psparam.psgap; emitlength = psgap/pslength;
    pstxelen = emitvec(1) : psgap : emitvec(end) - psgap;
    ec = 1;
    for txcount = 1:emitlength:length(emitnvec)-1
        ps = zeros(zlength , xlength);
        for framecount = txcount : txcount + emitlength-1
            pc = pstemp{framecount};
            ps = ps + pc; 
        end
        totalps{ec} = ps;
        ec = ec + 1
    end
    varargout{1} = totalps;
    end
end