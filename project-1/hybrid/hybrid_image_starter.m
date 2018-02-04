close all; % closes all figures

% read images and convert to single format
im1 = im2single(imread('./anakin_01c.jpg'));
im2 = im2single(imread('./anakin_01b.jpg'));
% im1 = rgb2gray(im1); % convert to grayscale
% im2 = rgb2gray(im2);

% use this if you want to align the two images (e.g., by the eyes) and crop
% them to be of same size
[im1, im2] = align_images(im1, im2);

% uncomment this when debugging hybridImage so that you don't have to keep aligning
% keyboard; 

%% Choose the cutoff frequencies and compute the hybrid image (you supply
%% this code)
close all;

R1 = im1(:,:,1);
G1 = im1(:,:,2);
B1 = im1(:,:,3);

R2 = im2(:,:,1);
G2 = im2(:,:,2);
B2 = im2(:,:,3);

s1 = 6;
r1 = 15;
s2 = 5;
r2 = 17;

R12 = hybridImage(R1, R2, 3, 9, 3, 9);
G12 = hybridImage(G1, G2, 3, 9, 3, 9);
B12 = hybridImage(B1, B2, 3, 9, 3, 9);

H = cat(3,R12,G12,B12);

H = rgb2gray(H);
imshow(H,[]);


%% Crop resulting image (optional)
im12=H;
figure(1), hold off, imagesc(im12), axis image, colormap gray
disp('input crop points');
[x, y] = ginput(2);  x = round(x); y = round(y);
im12 = im12(min(y):max(y), min(x):max(x), :);
figure(1), hold off, imagesc(im12), axis image, colormap gray