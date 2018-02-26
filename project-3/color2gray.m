function [im_out] = color2gray(imin)
   
% imin = imresize(im2double(imread('./samples/testing_gray.png')), 0.5, 'bilinear');
[imh,imw,~] = size(imin);
im2var = zeros(imh,imw);
im2var(1:imh*imw) = 1:imh* imw;
numeqs = 2*imh*imw + 1;

combo = zeros(imh,imw);

R = imin(:,:,1);
G = imin(:,:,2);
B = imin(:,:,3);

gray_val = imin(1,1,1);

A = sparse([], [], [], numeqs, imh*imw);
b = zeros(numeqs,1);

e=1;
for y=1:imh
    
    for x=1:imw-1
        disp(y);
        disp(x);
        A(e, im2var(y,x+1))=1; 
        A(e, im2var(y,x))=-1;
        
        
        b(e) = mymax( R(y,x+1)-R(y,x),  G(y,x+1)-G(y,x),  B(y,x+1)-B(y,x));
        
        e = e+1;
    end
end

for y=1:imh-1
    for x=1:imw
        A(e, im2var(y+1,x))=1; 
        A(e, im2var(y,x))=-1; 
        b(e) = mymax( R(y+1,x)-R(y,x),  G(y+1,x)-G(y,x),  B(y+1,x)-B(y,x));
%         b(e) = combo(y+1,x)-combo(y,x);
        e = e+1;
    end
end

A(e, im2var(1,1)) = 1;
b(e) = gray_val;

disp("Solving");
v = A\b;
disp("Solved");
im_out  =  reshape(v, [imh imw]);
% figure; imshow(im_out);



end