function[p_mask] = mk_linear_array(xpos, N_elements, pitch, Nx, Ny, dy)
    p_mask = zeros(Nx, Ny);
    ypos = ((1:N_elements)-(N_elements+1)/2)*pitch/dy+Ny/2;
    yvec = round(ypos);
    if (min(yvec)<1)
        error("grid too small for array!");
    end
    p_mask(xpos, yvec) = 1;
    [testx, testy] = size(p_mask);
    if(testx>Nx) || (testy>Ny)
        error("grid too small for array!");
    end
end