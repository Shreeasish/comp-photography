%% Load image
close all; % closes all figures

% read images and convert to single format
rgb = imread('low-contrast-3.jpg');
orig  = rgb;
R = rgb(:,:,1);
G = rgb(:,:,2);
B = rgb(:,:,3);

R = enhance(R);
G = enhance(G);
B = enhance(B);

rgb = cat(3,R,G,B);
%%

imshowpair(orig,uint8(rgb),'montage');


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
        M2(i,j) = floor(255*pn);
    end
end

end
