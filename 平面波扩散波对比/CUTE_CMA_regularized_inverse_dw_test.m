clear;
filename = './ps_info/ps_info_r3mm_is1450_bg1540_depth15.mat';
load(filename);

psvec = vectorize_ps(totalps , totalmask);
Mmat = calc_dw_CMA_mmat(totalmask , ROIparam , psparam , arrayparam);
phisize = [length(ROIparam.zvec) length(ROIparam.xvec)];

Dx = toeplitz([-1 1 zeros(1 ,phisize(2) -2)] , [-1 zeros(1 , phisize(2) -1)]);
Dx = kron(Dx(: , 1 : phisize(2) - 1).' , eye(phisize(1)));
Dz = toeplitz([-1 zeros(1 , phisize(1) -1)], [-1 1 zeros(1 , phisize(1) -2)]);
Dz = kron(eye(phisize(2)) , Dz(1 : phisize(1) - 1 , :));

gammax = 400e-6;
gammaz = 400e-6;


tic
Minv = (Mmat.' * Mmat + gammax * (Dx.' * Dx) + gammaz * (Dz.' * Dz)) \ Mmat.';
toc

tic
chat = Minv*psvec;
toc

ds=reshape(chat,phisize);
figure;
c=1540;
Cmap0 = 1./(1./c-ds);
imagesc(ROIparam.xvec*1e3, ROIparam.zvec*1e3, Cmap0);
colorbar
