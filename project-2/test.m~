I = im2double(imread('./samples/bricks_small.jpg'));
I = rgb2gray(I);

patchsize = 65;
numpatches = 4;
outputsize = patchsize*numpatches;
overlapsize = 16;

[s_y, s_x] = size(I);

synth = zeros(outputsize, outputsize);

rng(0,'twister')

rx = randi([1 s_x-patchsize],1,1);
ry = randi([1 s_y-patchsize],1,1);

patch = I(ry:ry+patchsize-1, rx:rx+patchsize-1,1);

mask  = cut(patch);