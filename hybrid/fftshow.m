function [F] = fftshow(img)
F = fftshift(img);
F = abs(F);
F = log(F+1);
imagesc(F);
end