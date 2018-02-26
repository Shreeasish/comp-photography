function [im_out] = toy_reconstruct(imin)
    [imh,imw,~] = size(imin);
    im2var = zeros(imh,imw);
    im2var(1:imh*imw) = 1:imh* imw;
    numeqs = 2*imh*imw + 1;

    A = sparse([], [], [], numeqs, imh*imw);
    b = zeros(numeqs,1);

    e=1;
    for y=1:imh
        for x=1:imw-1
            A(e, im2var(y,x+1))=1; 
            A(e, im2var(y,x))=-1; 
            b(e) = imin(y,x+1)-imin(y,x);
            e = e+1;
        end
    end

    for y=1:imh-1
        for x=1:imw
            A(e, im2var(y+1,x))=1; 
            A(e, im2var(y,x))=-1; 
            b(e) = imin(y+1,x)-imin(y,x);
            e = e+1;
        end
    end

    A(e, im2var(1,1)) = 1;
    b(e) = imin(1,1);

    v = A\b;
    im_out  =  reshape(v, [imh imw]);
    imshow(im_out);

end

