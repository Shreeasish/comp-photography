close all;

rgb = im2single(imread('mountain.jpg'));

lab = rgb2lab(rgb);
orig_rgb = rgb;

l = lab(:,:,1);
a = lab(:,:,2);
b = lab(:,:,3);

a_max = max(a(:));
a = a + a_max/2;

lab = cat(3,l,a,b);
out = lab2rgb(lab);
figure; imshow(out);
figure; imshow(orig_rgb);

imwrite(out,'./mountain_more_red.jpg');