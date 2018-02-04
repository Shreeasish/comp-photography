%% Load image
close all; % closes all figures

% read images and convert to single format
rgb = imread('low-contrast-4.jpg');
orig = rgb;

[h,s,v] = rgb2hsv(rgb);
v = floor(v*256);
v = enhance(v);
v = v/256;

hsv = cat(3,h,s,v);
rgb = hsv2rgb(hsv);

imshowpair(orig,rgb,'montage');


%%

function [M2] = enhance(M)

freq = zeros(256,1);

numPixels = size(M,1)*size(M,2);

Pn = zeros(256,1);

for i=1:size(M,1)
    for j=1:size(M,2)
        value=uint8(M(i,j));
        freq(value+1)=freq(value+1)+1;
        Pn(value+1) = freq(value+1)/numPixels;
    end
end


M2 = zeros(size(M,1),size(M,2));

for i=1:size(M,1)
    for j=1:size(M,2)
        pn = sum(Pn(1: uint8(M(i,j))+1 ));
        M2(i,j) = floor(255*pn);
    end
end

end
