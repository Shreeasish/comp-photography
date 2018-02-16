I = im2double(imread('./samples/fabric-small.jpg'));
% I = rgb2gray(I);

[h, s, I] = rgb2hsv(I);

patchsize = 65;
numpatches = 4;
% outputsize = 2;
% numpatches = outputsize/patchsize;
outputsize = patchsize*numpatches;
overlapsize = 16;

[s_y, s_x] = size(I);

synth = zeros(outputsize, outputsize);
synthh = zeros(outputsize, outputsize);
synths = zeros(outputsize, outputsize);

rng(0, 'twister')

rx = randi([1 s_x-patchsize],1,1);
ry = randi([1 s_y-patchsize],1,1);

patch = I(ry:ry+patchsize-1, rx:rx+patchsize-1,1);
spatch = s(ry:ry+patchsize-1, rx:rx+patchsize-1,1);
hpatch = h(ry:ry+patchsize-1, rx:rx+patchsize-1,1);

synth(1:patchsize,1:patchsize) = patch(:,:);
synths(1:patchsize,1:patchsize) = spatch(:,:);
synthh(1:patchsize,1:patchsize) = hpatch(:,:);

T = zeros(patchsize,patchsize);
M = zeros(patchsize,patchsize);
M(:,1:overlapsize) = 1;

for i=1:numpatches
    T(:,1:overlapsize) = synth(1 : patchsize, 1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize );
    ssd_matrix = ssd_patch(T,M,I);
    [patch,spatch,hpatch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x, s, h);
    synth( 1 : patchsize, 1 + i*patchsize - (i-1)*overlapsize - overlapsize/2: (i+1)*patchsize - i*overlapsize )  = patch(:, overlapsize/2 + 1 :end);
    synths( 1 : patchsize, 1 + i*patchsize - (i-1)*overlapsize - overlapsize/2: (i+1)*patchsize - i*overlapsize )  = spatch(:, overlapsize/2 + 1 :end);
    synthh( 1 : patchsize, 1 + i*patchsize - (i-1)*overlapsize - overlapsize/2: (i+1)*patchsize - i*overlapsize )  = hpatch(:, overlapsize/2 + 1 :end);
end

%%1 : patchsize,
for i=1:numpatches
    for j=1:numpatches
        T = zeros(patchsize,patchsize);
        M = zeros(patchsize,patchsize);
        if j == 1
           T(1:overlapsize,:) = synth( 1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize , 1 : patchsize);
           M(1:overlapsize,:) = 1;
           ssd_matrix = ssd_patch(T,M,I);
           [patch, spatch, hpatch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x, s, h);
           synth( 1 + i*patchsize - (i-1) * overlapsize - overlapsize/2 :  (i+1)*patchsize - i*overlapsize , 1: patchsize) = patch(overlapsize/2 + 1 : end, :);
           synths( 1 + i*patchsize - (i-1) * overlapsize - overlapsize/2 :  (i+1)*patchsize - i*overlapsize , 1: patchsize) = spatch(overlapsize/2 + 1 : end, :);
           synthh( 1 + i*patchsize - (i-1) * overlapsize - overlapsize/2 :  (i+1)*patchsize - i*overlapsize , 1: patchsize) = hpatch(overlapsize/2 + 1 : end, :);
        end
            T(1:overlapsize,1:patchsize) ...
              = synth( 1 +  i*patchsize - i*overlapsize : i*patchsize - (i-1)*overlapsize , 1 + j*patchsize - j*overlapsize : (j+1)*patchsize -j*overlapsize);
            T(1:patchsize , 1 : overlapsize) ...
                = synth( 1 + i*patchsize - i*overlapsize : (i+1)*patchsize -i*overlapsize , 1 +  j*patchsize - j*overlapsize : j*patchsize - (j-1)*overlapsize);
            M(1:overlapsize,:) = 1;
            M(:,1:overlapsize) = 1;
            [patch, spatch, hpatch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x, s, h);
            synth( 1 + i*patchsize - (i-1) * overlapsize - overlapsize/2 :  (i+1)*patchsize - i*overlapsize , ...
                1 + j*patchsize - (j-1)*overlapsize - overlapsize/2: (j+1)*patchsize - j*overlapsize) = patch(overlapsize/2 + 1 :end, overlapsize/2 + 1 :end);
            
            synths( 1 + i*patchsize - (i-1) * overlapsize - overlapsize/2 :  (i+1)*patchsize - i*overlapsize , ...
                1 + j*patchsize - (j-1)*overlapsize - overlapsize/2: (j+1)*patchsize - j*overlapsize) = spatch(overlapsize/2 + 1 :end, overlapsize/2 + 1 :end);
            
            synthh( 1 + i*patchsize - (i-1) * overlapsize - overlapsize/2 :  (i+1)*patchsize - i*overlapsize , ...
                1 + j*patchsize - (j-1)*overlapsize - overlapsize/2: (j+1)*patchsize - j*overlapsize) = hpatch(overlapsize/2 + 1 :end, overlapsize/2 + 1 :end);
    end
end

synth = cat(3,synthh,synths,synth);
synth = hsv2rgb(synth);

imshow(synth,[]);

imwrite(synth,'./results/fabric-overlap.jpg');


function [patch, spatch, hpatch] = choose_sample(ssd_matrix, I, patchsize, s_y, s_x, s, h)
    pby2 = floor(patchsize/2);
    [R,C] = ndgrid( 1 : size(ssd_matrix,1), 1:size(ssd_matrix,2) );
    [sorted_ssd, idx] = sort(ssd_matrix(:));
    r = randi([1 400],1,1);
    while not(( R(idx(r)) - pby2 > 0 ) && ( R(idx(r)) + pby2 < s_y) && ( C(idx(r))- pby2> 0 ) && (C(idx(r)) + pby2 < s_x ))
        r = randi([1 400],1,1);
        disp(1);
    end
    patch = I( R(idx(r)) - pby2 : R(idx(r)) + pby2 , C(idx(r))- pby2 : C(idx(r)) + pby2 );
    spatch = s( R(idx(r)) - pby2 : R(idx(r)) + pby2 , C(idx(r))- pby2 : C(idx(r)) + pby2 );
    hpatch = h( R(idx(r)) - pby2 : R(idx(r)) + pby2 , C(idx(r))- pby2 : C(idx(r)) + pby2 );
% patch =0;
end

function[ssd] = ssd_patch(T,M,I)
    ssd = imfilter(I.^2, M) -2*imfilter(I, M.*T) + sum(sum((M.*T).^2));
end