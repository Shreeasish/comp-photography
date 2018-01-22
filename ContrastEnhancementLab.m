%% Load image
close all; % closes all figures

% read images and convert to single format
rgb = imread('low-contrast.jpg');
orig  = rgb;
lab = rgb2lab(rgb);
l = lab(:,:,1);
a = lab(:,:,2);
b = lab(:,:,3);
%%

maxL = max(l(:));
maxA = max(a(:)); 
maxB = max(b(:));

minL = min(l(:));
minA = min(a(:)); 
minB = min(b(:));

a = (((a - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin



l = enhance(l);
a = enhance(a);
b = enhance(b);

lab = cat(3,l,a,b);
%%
rgb = lab2rgb(lab);
imshowpair(orig,rgb,'montage');

%%

function [M2] = enhance(M)


freq = zeros(256,1);

numPixels = size(M,1)*size(M,2);

Pn = zeros(256,1);

for i=1:size(M,1)
    for j=1:size(M,2)
        value=round(M(i,j));
        freq(value+1)=freq(value+1)+1;
        Pn(value+1) = freq(value+1)/numPixels;
    end
end


M2 = zeros(size(M,1),size(M,2));

for i=1:size(M,1)
    for j=1:size(M,2)
        pn = sum(Pn(1: round(M(i,j))+1 ));
        M2(i,j) = floor(128*pn);
    end
end

end
