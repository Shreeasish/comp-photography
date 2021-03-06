input_texture = imread('./samples/fabric-small.jpg');
patchsize = 65;
OUTPUTSIZE = 261;
num_patches = 4;

[s_y, s_x, s_z] = size(input_texture);
synth = zeros(OUTPUTSIZE,OUTPUTSIZE,3);

for i=1:num_patches
    for j=1:num_patches
        rx = randi([1 s_x-patchsize],1,1);
        ry = randi([1 s_y-patchsize],1,1);
        patch1 = input_texture(ry:ry+patchsize-1, rx:rx+patchsize-1,1);
        patch2 = input_texture(ry:ry+patchsize-1, rx:rx+patchsize-1,2);
        patch3 = input_texture(ry:ry+patchsize-1, rx:rx+patchsize-1,3);
        X = (i-1)*patchsize + 1;
        pi = 1;
        for X=X:X+patchsize-1
            pj = 1;
            Y = (j-1)*patchsize + 1;
            for Y=Y:Y+patchsize-1
                synth(X,Y,1) = patch1(pi,pj);
                synth(X,Y,2) = patch2(pi,pj);
                synth(X,Y,3) = patch3(pi,pj);
                pj = pj + 1;
            end
            pi = pi + 1;
        end
        
    end
end

imshow(uint8(synth));
imwrite(uint8(synth),'./results/fabric-basic.jpg');


