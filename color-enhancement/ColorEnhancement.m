close all;
rgb = imread('./color-enhancement/bad_color-2.jpg');

orig_rgb = rgb;

lab = rgb2lab(rgb);

l = lab(:,:,1);
a = lab(:,:,2);
b = lab(:,:,3);

a = a*2;
b = b*2;

lab = cat(3,l,a,b);
out = lab2rgb(lab);
figure; imshowpair(out,orig_rgb,'montage');
imwrite(out,'./color-enhancement/color_corrected_2.jpg');