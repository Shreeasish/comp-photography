close all;

rgb = im2single(imread('mountain.jpg'));

lab = rgb2lab(rgb);
orig_rgb = rgb;

l = lab(:,:,1);
a = lab(:,:,2);
b = lab(:,:,3);

b_max = max(b(:));
b = b - b_max/2;

lab = cat(3,l,a,b);
out = lab2rgb(lab);
figure; imshow(out);
figure; imshow(orig_rgb);
imwrite(out,'./mountain_less_yellow.jpg');