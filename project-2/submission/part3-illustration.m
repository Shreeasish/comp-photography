clear;
close all;

I = im2double(imread('./samples/fabric-small.jpg'));
% I = rgb2gray(I);

[h, s, I] = rgb2hsv(I);

patchsize = 65;
numpatches = 5;
overlapsize = 16;
outputsize = patchsize*numpatches - (numpatches-1)*overlapsize;

[s_y, s_x] = size(I);

synth = zeros(outputsize, outputsize);
synthh = zeros(outputsize, outputsize);
synths = zeros(outputsize, outputsize);

rng(0,'twister')

rx = randi([1 s_x-patchsize],1,1);
ry = randi([1 s_y-patchsize],1,1);

hpatch = h(ry:ry+patchsize-1, rx:rx+patchsize-1,1);
spatch = s(ry:ry+patchsize-1, rx:rx+patchsize-1,1);
patch = I(ry:ry+patchsize-1, rx:rx+patchsize-1,1);


synthh(1:patchsize,1:patchsize) = hpatch(:,:);
synths(1:patchsize,1:patchsize) = spatch(:,:);
synth(1:patchsize,1:patchsize) = patch(:,:);

T = zeros(patchsize,patchsize);
M = zeros(patchsize,patchsize);
M(:,1:overlapsize) = 1;

for i=1:2
    T(:,1:overlapsize) = synth(1 : patchsize, 1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize );
    ssd_matrix = ssd_patch(T,M,I);
    
    [hpatch, spatch, patch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x,h,s);
    
    im2 = cat(3,hpatch,spatch,patch);
    im2 = hsv2rgb(im2);
    imwrite(im2,'./results/patch2.jpg')
    
    
    im1 = cat(3,synthh( 1: patchsize, 1 + (i-1)*patchsize - (i-1)*overlapsize : (i)*patchsize - (i-1)*overlapsize ), ...
        synths( 1: patchsize, 1 + (i-1)*patchsize - (i-1)*overlapsize : (i)*patchsize - (i-1)*overlapsize ), ...
        synth( 1: patchsize, 1 + (i-1)*patchsize - (i-1)*overlapsize : (i)*patchsize - (i-1)*overlapsize ));
    im1= hsv2rgb(im1);
    imwrite(im1,'./results/patch1.jpg')
    
    overlap1 = synth(1 : patchsize, 1 + i*patchsize -i*overlapsize : i*patchsize - (i-1)*overlapsize);
    overlap2 = patch(1 : patchsize, 1:overlapsize);
    
    
    
    bndcost = (overlap1-overlap2).^2;
    [mask, bestpath] = cut(transpose(bndcost));
    mask = transpose(mask);
    
    hpatch(:,1:overlapsize) = double(mask).*hpatch(:,1:overlapsize);
    spatch(:,1:overlapsize) = double(mask).*spatch(:,1:overlapsize);
    patch(:,1:overlapsize) = double(mask).*patch(:,1:overlapsize);
    
    im1 = cat(3,synthh( 1: patchsize, 1 + (i-1)*patchsize - (i-1)*overlapsize : (i)*patchsize - (i-1)*overlapsize ), ...
        synths( 1: patchsize, 1 + (i-1)*patchsize - (i-1)*overlapsize : (i)*patchsize - (i-1)*overlapsize ), ...
        synth( 1: patchsize, 1 + (i-1)*patchsize - (i-1)*overlapsize : (i)*patchsize - (i-1)*overlapsize ));
    im1= hsv2rgb(im1);
    
    imwrite(im1,'./results/patch1.jpg')
    
    
    
    
    
    normbndcost = bndcost - min(bndcost(:));
    normbndcost = normbndcost ./ max(normbndcost(:));
    imshow(normbndcost); axis image; hold on; plot(bestpath, 1:length(bestpath));
    
    synthh(1:patchsize,1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize) ...
        = ~mask.*synthh(1:patchsize,1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize );
    synths(1:patchsize,1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize) ...
        = ~mask.*synths(1:patchsize,1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize );
    synth(1:patchsize,1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize) ...
        = ~mask.*synth(1:patchsize,1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize );
        
    synthh( 1: patchsize, 1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize ) ...
        = synthh( 1: patchsize, 1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize ) + hpatch;
    synths( 1: patchsize, 1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize ) ...
        = synths( 1: patchsize, 1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize ) + spatch;
    synth( 1: patchsize, 1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize ) ...
        = synth( 1: patchsize, 1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize ) + patch;
end

synth = cat(3,synthh,synths,synth);
synth = hsv2rgb(synth);

%%

for i=1:numpatches -1 
    for j=1:numpatches - 1
        T = zeros(patchsize,patchsize);
        M = zeros(patchsize,patchsize);
        if j == 1
            T(1:overlapsize,:) = synth( 1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize , 1 : patchsize);
            M(1:overlapsize,:) = 1;
            ssd_matrix = ssd_patch(T,M,I);
            
            [hpatch, spatch, patch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x,h,s);

            overlap1 = synth( 1 + i*patchsize -i*overlapsize : i*patchsize - (i-1)*overlapsize, 1:patchsize);

            overlap2 = patch( 1:overlapsize, 1 : patchsize);

            bndcost = (overlap1-overlap2).^2;
            mask = cut(bndcost);
            
            synthh(1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize, 1:patchsize) ...
                = ~mask.*synthh(1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize, 1:patchsize);
            synths(1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize, 1:patchsize) ...
                = ~mask.*synths(1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize, 1:patchsize);
            synth(1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize, 1:patchsize) ...
                = ~mask.*synth(1 + i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize, 1:patchsize);
    
            hpatch(1:overlapsize,:) = double(mask).*hpatch(1:overlapsize,:);
            spatch(1:overlapsize,:) = double(mask).*spatch(1:overlapsize,:);
            patch(1:overlapsize,:) = double(mask).*patch(1:overlapsize,:);
            
            synthh( 1 + i*patchsize - i*overlapsize :  (i+1)*patchsize - i*overlapsize , ...
                1: patchsize) = synthh( 1 + i*patchsize - i*overlapsize :  (i+1)*patchsize - i*overlapsize , 1: patchsize) + hpatch;
            synths( 1 + i*patchsize - i*overlapsize :  (i+1)*patchsize - i*overlapsize , ...
                1: patchsize) = synths( 1 + i*patchsize - i*overlapsize :  (i+1)*patchsize - i*overlapsize , 1: patchsize) + spatch;
            synth( 1 + i*patchsize - i*overlapsize :  (i+1)*patchsize - i*overlapsize , ...
                1: patchsize) = synth( 1 + i*patchsize - i*overlapsize :  (i+1)*patchsize - i*overlapsize , 1: patchsize) + patch;

        end
            T(1:overlapsize,1:patchsize) ...
              = synth( 1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize , 1 + j*patchsize - j*overlapsize : (j+1)*patchsize -j*overlapsize);
            T(1:patchsize , 1 : overlapsize) ...
                = synth( 1 + i*patchsize - i*overlapsize : (i+1)*patchsize -i*overlapsize , 1 +  j*patchsize - j*overlapsize : j*patchsize - (j-1)*overlapsize);
            M(1:overlapsize,:) = 1;
            M(:,1:overlapsize) = 1;
            
            ssd_matrix = ssd_patch(T,M,I);
            
            [hpatch, spatch, patch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x,h,s);
            
            overlap_h1 = ...
              synth( 1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize , 1 + j*patchsize - j*overlapsize : (j+1)*patchsize -j*overlapsize);
            
            overlap_h2 = ...
                patch(1:overlapsize,1:patchsize);

            bndcost = (overlap_h1-overlap_h2).^2;
            mask1 = cut(bndcost);
            
            
            overlap_v1 = ...
                synth( 1 + i*patchsize - i*overlapsize : (i+1)*patchsize -i*overlapsize , 1 +  j*patchsize - j*overlapsize : j*patchsize - (j-1)*overlapsize);
            
            overlap_v2 = ...
                patch(1:patchsize,1:overlapsize);

            bndcost = (overlap_v1-overlap_v2).^2;
            mask2 = transpose(cut(transpose(bndcost)));
            
            
            mask1(1:overlapsize,1:overlapsize) = mask1(1:overlapsize,1:overlapsize) & mask2(1:overlapsize,1:overlapsize);
            mask2(1:overlapsize,1:overlapsize) = mask1(1:overlapsize,1:overlapsize) & mask2(1:overlapsize,1:overlapsize);
     
            
            hpatch(1:overlapsize,:) = mask1.*hpatch(1:overlapsize,:);
            hpatch(:,1:overlapsize) = mask2.*hpatch(:,1:overlapsize);
            
            spatch(1:overlapsize,:) = mask1.*spatch(1:overlapsize,:);
            spatch(:,1:overlapsize) = mask2.*spatch(:,1:overlapsize);
            
            patch(1:overlapsize,:) = mask1.*patch(1:overlapsize,:);
            patch(:,1:overlapsize) = mask2.*patch(:,1:overlapsize);
            
            
            synthh(1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize , ...
                1 + j*patchsize - j*overlapsize : (j+1)*patchsize - j*overlapsize) ...
                = ~mask1.*synthh(1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize, ...
                1 + j*patchsize - j*overlapsize : (j+1)*patchsize - j*overlapsize);
            
            synths(1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize , ...
                1 + j*patchsize - j*overlapsize : (j+1)*patchsize - j*overlapsize) ...
                = ~mask1.*synths(1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize, ...
                1 + j*patchsize - j*overlapsize : (j+1)*patchsize - j*overlapsize);
            
            synth(1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize , ...
                1 + j*patchsize - j*overlapsize : (j+1)*patchsize - j*overlapsize) ...
                = ~mask1.*synth(1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize, ...
                1 + j*patchsize - j*overlapsize : (j+1)*patchsize - j*overlapsize);
            
            
            
            synthh(1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize, ...
                1 + j*patchsize - j*overlapsize : j*patchsize - (j-1)*overlapsize ) ...
                = ~mask2.*synthh(1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize, ...
                1 + j*patchsize - j*overlapsize : j*patchsize - (j-1)*overlapsize );
            
            synths(1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize, ...
                1 + j*patchsize - j*overlapsize : j*patchsize - (j-1)*overlapsize ) ...
                = ~mask2.*synths(1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize, ...
                1 + j*patchsize - j*overlapsize : j*patchsize - (j-1)*overlapsize );
            
            synth(1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize, ...
                1 + j*patchsize - j*overlapsize : j*patchsize - (j-1)*overlapsize ) ...
                = ~mask2.*synth(1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize, ...
                1 + j*patchsize - j*overlapsize : j*patchsize - (j-1)*overlapsize );
            
            
            
            synthh( 1 + i*patchsize - i * overlapsize:  (i+1)*patchsize - i*overlapsize , ...
                1 + j*patchsize - j*overlapsize: (j+1)*patchsize - j*overlapsize) ...
                = synthh( 1 + i*patchsize - i * overlapsize:  (i+1)*patchsize - i*overlapsize , ...
                1 + j*patchsize - j*overlapsize: (j+1)*patchsize - j*overlapsize) ... 
                + hpatch;
            
            synths( 1 + i*patchsize - i * overlapsize:  (i+1)*patchsize - i*overlapsize , ...
                1 + j*patchsize - j*overlapsize: (j+1)*patchsize - j*overlapsize) ...
                = synths( 1 + i*patchsize - i * overlapsize:  (i+1)*patchsize - i*overlapsize , ...
                1 + j*patchsize - j*overlapsize: (j+1)*patchsize - j*overlapsize) ... 
                + spatch;
            
            synth( 1 + i*patchsize - i * overlapsize:  (i+1)*patchsize - i*overlapsize , ...
                1 + j*patchsize - j*overlapsize: (j+1)*patchsize - j*overlapsize) ...
                = synth( 1 + i*patchsize - i * overlapsize:  (i+1)*patchsize - i*overlapsize , ...
                1 + j*patchsize - j*overlapsize: (j+1)*patchsize - j*overlapsize) ... 
                + patch;
    end
end

synth = cat(3,synthh,synths,synth);
synth = hsv2rgb(synth);
imshow(synth,[]);

imwrite(synth,'./results/fabric-seam.jpg');


function [hpatch, spatch, patch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x, h, s)
    pby2 = floor(patchsize/2);
    [R,C] = ndgrid( 1 : size(ssd_matrix,1), 1:size(ssd_matrix,2) );
    [sorted_ssd, idx] = sort(ssd_matrix(:));
    r = randi([1 100],1,1);
    while not(( R(idx(r)) - pby2 > 0 ) && ( R(idx(r)) + pby2 < s_y) && ( C(idx(r))- pby2> 0 ) && (C(idx(r)) + pby2 < s_x ))
        r = randi([1 100],1,1);
    end
    
    disp("Selected");
    
    patch = I( R(idx(r)) - pby2 : R(idx(r)) + pby2 , C(idx(r))- pby2 : C(idx(r)) + pby2 );
    spatch = s( R(idx(r)) - pby2 : R(idx(r)) + pby2 , C(idx(r))- pby2 : C(idx(r)) + pby2 );
    hpatch = h( R(idx(r)) - pby2 : R(idx(r)) + pby2 , C(idx(r))- pby2 : C(idx(r)) + pby2 );
% patch =0;
end


function[ssd] = ssd_patch(T,M,I)
    ssd = imfilter(I.^2, M) -2*imfilter(I, M.*T) + sum(sum((M.*T).^2));
end