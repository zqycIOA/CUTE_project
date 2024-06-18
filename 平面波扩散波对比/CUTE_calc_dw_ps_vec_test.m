clear;
% close all;
fs = 50e6;
sigmaxvec = (-20 : 0.5 : 20) * 1e-3; sigmazvec = (5 : 0.5 : 40) * 1e-3;
[sigmaxs , sigmazs] = meshgrid(sigmaxvec , sigmazvec);
pitch = 0.300e-3;
N_elements = 128;
axiscor = linspace(-0.5, 0.5, N_elements)*pitch*(N_elements-1);

axiscortx = (1 : N_elements) * pitch;
axiscortx = axiscortx - mean(axiscortx);
axiscorrx = axiscortx(1):pitch/2:axiscortx(end);
arrayparam.axiscortx = axiscortx;
arrayparam.axiscorrx = axiscorrx;

ROIparam.xvec = sigmaxvec;
ROIparam.zvec = sigmazvec;
ROIparam.c = 1540;

psparam.CMAvec = [-15 0 15] * pi / 180;
psparam.psstart = 48;
psparam.psend = 80;
psparam.psgap = 8;
psparam.pslength = 4;

linevec = 1:128;
lc = 1;
for linecount = linevec
    filename = ['./rf_signal/dw_r3mm_x0z15_is1450_bg1540_255rx/line_' num2str(linecount) '.mat'];
    load(filename);
    rf_frame.emitele = linecount;
    totalframes(lc) = rf_frame;
    clear rf_frame; 
    lc = lc + 1;
end

[totalps , totalmask] = calc_phase_shift_vec_dw_CMA(totalframes , psparam , ROIparam , arrayparam);
% totalps = calc_phase_shift_vec_dw_CMA(totalframes , psparam , ROIparam , arrayparam);

inclusion_map = sqrt(sigmaxs.^2 + (sigmazs - 15e-3).^2) < 3e-3;
ROIparam.Cmap0 = 1540 * (1 - inclusion_map) + 1450 * inclusion_map;
psref = calc_dw_ps_prediction(ROIparam , psparam , arrayparam);

figure;
imagesc(sigmaxvec * 1e3 , sigmazvec * 1e3 , totalps{1,2});
caxis([-2e-7 2e-7]); colorbar
axis([-20 20 0 40])

figure;
imagesc(sigmaxvec * 1e3 , sigmazvec * 1e3 , psref{1,2}.*totalmask{1,2});
caxis([-2e-7 2e-7]); colorbar
axis([-20 20 0 40])

filename = './ps_info/ps_info_r3mm_is1450_bg1540_depth15_tbiased.mat';
save(filename , 'totalmask' , 'totalps' , 'ROIparam' , 'psparam' , 'arrayparam');