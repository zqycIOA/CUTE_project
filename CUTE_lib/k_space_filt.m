function[I_filt] = k_space_filt(I , xvec , zvec , angle)
    Ikxkz = fftshift(fft2(I));
    dx = xvec(2) - xvec(1);
    dz = zvec(2) - zvec(1);
    [nz , nx] = size(I);
    kxvec = linspace(-0.5/dx , 0.5/dx , nx);
    kzvec = linspace(-0.5/dz , 0.5/dz , nz);
    [kxs , kzs] = meshgrid(kxvec , kzvec);
    kfilt = exp(-abs(atan(kxs./kzs) - angle).^2/(10 * pi / 180).^2);
    kfilt(isnan(kfilt)) = 1;
    Ikxkz_filt = ifftshift(kfilt .* Ikxkz);
    I_filt = ifft2(Ikxkz_filt);
end