clear;
rng(0,'twister')

alpha = 0.1;

I = im2double(imread('./samples/frost.jpg'));
[h, s, I] = rgb2hsv(I);
[s_y, s_x] = size(I);

original_target  = im2double(imread('./samples/low-contrast-4.jpg'));


graytarget = rgb2gray(original_target);
original_target_size = size(graytarget);

i = imgaussfilt(I,3);
graytarget = imgaussfilt(graytarget,3);

patchsize = 13;
overlapsize = 4;

numpatches_x = floor(original_target_size(2)/(patchsize));
numpatches_y = floor(original_target_size(1)/(patchsize));

synthh = zeros(original_target_size(1),original_target_size(2));
synths = zeros(original_target_size(1),original_target_size(2));
synth = zeros(original_target_size(1),original_target_size(2));


M2 = ones(patchsize, patchsize);
T2 = graytarget(1:patchsize, 1:patchsize);

ssd2 = ssd_patch(T2,M2,I)*(1-alpha);

[hpatch, spatch, patch] = choose_sample(ssd2, I, patchsize, s_y, s_x, h,s);

synthh(1:patchsize, 1:patchsize) = hpatch;
synths(1:patchsize, 1:patchsize) = spatch;
synth(1:patchsize, 1:patchsize) = patch;

T = zeros(patchsize,patchsize);
M = zeros(patchsize,patchsize);
M(:,1:overlapsize) = 1;


for i=1:numpatches_x + 8 
    T(:,1:overlapsize) = synth(1 : patchsize, 1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize );
    T2 = graytarget(1:patchsize, 1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize );
    
    ssd_matrix = ssd_patch(T,M,I)*alpha;
    ssd2 = ssd_patch(T2,M2,I)*(1-alpha);
    ssd_matrix = ssd_matrix + ssd2;
    
    [hpatch, spatch, patch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x,h,s);
    
    overlap1 = synth(1 : patchsize, 1 + i*patchsize -i*overlapsize : i*patchsize - (i-1)*overlapsize);
    
    overlap2 = patch(1 : patchsize, 1:overlapsize);
    
    bndcost = (overlap1-overlap2).^2;
    mask = transpose(cut(transpose(bndcost)));
    
    hpatch(:,1:overlapsize) = double(mask).*hpatch(:,1:overlapsize);
    spatch(:,1:overlapsize) = double(mask).*spatch(:,1:overlapsize);
    patch(:,1:overlapsize) = double(mask).*patch(:,1:overlapsize);
    
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


for i=1:numpatches_y + 8
    for j=1:numpatches_x + 8
        disp(i);
        disp(j);
        T = zeros(patchsize,patchsize);
        M = zeros(patchsize,patchsize);
        if j == 1
            T(1:overlapsize,:) = synth( 1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize , 1 : patchsize);   
            M(1:overlapsize,:) = 1;
         
            T2 = graytarget( 1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize, 1:patchsize );
            
            ssd_matrix = ssd_patch(T,M,I)*alpha;
            ssd2 = ssd_patch(T2,M2,I)*(1-alpha);
            
            ssd_matrix = ssd_matrix + ssd2;
            
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
            
            T2 = graytarget(1 + i*patchsize - i*overlapsize : (i+1)*patchsize - i*overlapsize,...
                1 + j*patchsize -  j*overlapsize : (j+1)*patchsize - j*overlapsize );
            
            ssd_matrix = ssd_patch(T,M,I)*alpha;
            ssd2 = ssd_patch(T2,M2,I)*(1-alpha);
            
            ssd_matrix = ssd_matrix + ssd2;
            
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


function [hpatch, spatch, patch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x, h, s)
    pby2 = floor(patchsize/2);
    [R,C] = ndgrid( 1 : size(ssd_matrix,1), 1:size(ssd_matrix,2) );
    [sorted_ssd, idx] = sort(ssd_matrix(:));
    r = randi([1 1000],1,1);
    while not(( R(idx(r)) - pby2 > 0 ) && ( R(idx(r)) + pby2 < s_y) && ( C(idx(r))- pby2> 0 ) && (C(idx(r)) + pby2 < s_x ))
        r = randi([1 1200],1,1);
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