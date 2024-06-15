function[psvec] = vectorize_ps(totalps , totalmask)
    [nCMA , nEle] = size(totalps);
    pssize = 0;
    for cc = 1 : nCMA
        for ec = 1 : nEle
            pssize = pssize + sum(sum(totalmask{cc , ec}));
        end
    end
    psvec = zeros(pssize , 1);
    pc = 1;
    for cc = 1 : nCMA
        for ec = 1 : nEle
            pstemp = totalps{cc , ec};
            masktemp = totalmask{cc , ec};
            masksize = sum(sum(masktemp));
            pstemp_masked = pstemp(masktemp == 1); 
            psvec(pc : pc + masksize -1,1) = pstemp_masked;
            pc = pc + masksize;
        end
    end
end