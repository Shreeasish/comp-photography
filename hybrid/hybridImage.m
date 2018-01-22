function [hybrid] = hybridImage(im1,im2,s1,r1,s2,r2)

g = compute_kernel(s1,r1);
g2 = compute_kernel(s2,r2);

im1_low = conv2(im1,g,'same');
im2_low = conv2(im2,g2,'same');

F = fft2(im2);
F_low = fft2(im2_low);
F_high = F - F_low;

IM1_low = fft2(im1_low);



HYBRID = IM1_low + F_high;


% HYBRID = imnoise(HYBRID, 'gaussian', 0, v / 10);


hybrid = ifft2(HYBRID);



% figure; imshow(hybrid,[]);
end

