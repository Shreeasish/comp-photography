function [im_final] = mixedBlend(im_s, mask_s, im_background)
%     im_background = imresize(im2double(imread('./samples/im2.jpg')), 0.25, 'bilinear');
%     im_object = imresize(im2double(imread('./samples/penguin-chick.jpeg')), 0.25, 'bilinear');
% 
%     % get source region mask from the user
%     objmask = getMask(im_object);
%     % align im_s and mask_s with im_background
%     [im_s, mask_s] = alignSource(im_object, objmask, im_background);

    [j,k,~] = find(rgb2gray(im_s));
    j_beg = j(1);
    k_beg = k(1);
    j_end = j(end);
    k_end = k(end);

    imh = j_end - j_beg + 1;
    imw = k_end - k_beg + 1;

    im2var = zeros(imh,imw);
    im2var(1:imh*imw) = 1:imh* imw;
    numeqs = 4*imh*imw;

    v = {};
    
    alfa = 0.5;
    

    for c=1:3


        Y=1;
        e = 1;

        A = sparse([], [], [], numeqs, imh*imw);
        b = zeros(numeqs,1);

        for y=j_beg:j_end
            X=1;
            for x=k_beg:k_end

                if ~(mask_s(y,x))

                    A(e, im2var(Y,X)) = 1;
                    b(e) = im_background(y,x,c);
                    e=e+1;

                    A(e, im2var(Y,X)) = 1;
                    b(e) = im_background(y,x,c);
                    e=e+1;

                    A(e, im2var(Y,X)) = 1;
                    b(e) = im_background(y,x,c);
                    e=e+1;

                    A(e, im2var(Y,X)) = 1;
                    b(e) = im_background(y,x,c);
                    e=e+1;
                else


                    % Up (y-1)
                    if(mask_s(y-1,x))
                        A(e,im2var(Y,X)) = 1;
                        A(e,im2var(Y-1,X)) = -1;
                        
%                         b(e) = im_s(y,x,c) - im_s(y-1,x,c);

                        s = im_s(y,x,c) - im_s(y-1,x,c);
                        t = im_background(y,x,c) - im_background(y-1,x,c);
                        d = alfa*s + (1-alfa)*t;
                        
                        b(e) = d;
                    else
                        A(e,im2var(Y,X)) = 1;
    %                     A(e,im2var(Y-1,X)) = -1;
%                         b(e) = im_s(y,x,c) - im_s(y-1,x,c) + im_background(y-1,x,c);

                        s = im_s(y,x,c) - im_s(y-1,x,c);
                        t = im_background(y,x,c) - im_background(y-1,x,c);
                        d = alfa*s + (1-alfa)*t ;
                        
                        b(e) = d + im_background(y-1,x,c);
                    end
                    e=e+1;

                    % Down (y+1)
                    if(mask_s(y+1,x))
                        A(e,im2var(Y,X)) = 1;
                        A(e,im2var(Y+1,X)) = -1;
%                         b(e) = im_s(y,x,c) - im_s(y+1,x,c);
                        
                        s = im_s(y,x,c) - im_s(y+1,x,c);
                        t = im_background(y,x,c) - im_background(y+1,x,c);
                        d = alfa*s + (1-alfa)*t;
                        
                        b(e) = d;
                    else
                        A(e,im2var(Y,X)) = 1;
    %                     A(e,im2var(Y+1,X)) = -1;
%                         b(e) = im_s(y,x,c) - im_s(y+1,x,c) + im_background(y+1,x,c);
                        s = im_s(y,x,c) - im_s(y+1,x,c);
                        t = im_background(y,x,c) - im_background(y+1,x,c);
                        d = alfa*s + (1-alfa)*t;
                        
                        b(e) = d + im_background(y+1,x,c);
                    end
                    e=e+1;

                    % right (x+1)
                    if(mask_s(y,x+1))
                        A(e,im2var(Y,X)) = 1;
                        A(e,im2var(Y,X+1)) = -1;
                        
%                         b(e) = im_s(y,x,c) - im_s(y,x+1,c);
                        
                        s = im_s(y,x,c) - im_s(y+1,x,c);
                        t = im_background(y,x,c) - im_background(y,x+1,c);
                        d = alfa*s + (1-alfa)*t;
                        
                        b(e) = d;
                    else
                        A(e,im2var(Y,X)) = 1;
                        %  A(e,im2var(Y,X+1)) = -1;
                        %  b(e) = im_s(y,x,c) - im_s(y,x+1,c) + im_background(y,x+1,c);
                        
                        s = im_s(y,x,c) - im_s(y,x+1,c);
                        t = im_background(y,x,c) - im_background(y,x+1,c);
                        d = alfa*s + (1-alfa)*t;
                        
                        b(e) = d + im_background(y,x+1,c);
                    end
                    e=e+1;

                    % left (x-1)
                    if(mask_s(y,x-1))
                        A(e,im2var(Y,X)) = 1;
                        A(e,im2var(Y,X-1)) = -1;
                        
%                         b(e) = im_s(y,x,c) - im_s(y,x-1,c);
                        
                        s = im_s(y,x,c) - im_s(y,x-1,c);
                        t = im_background(y,x,c) - im_background(y,x-1,c);
                        d = alfa*s + (1-alfa)*t;
                        
                        b(e) = d;
                    else
                        A(e,im2var(Y,X)) = 1;
    %                     A(e,im2var(Y,X-1)) = -1;
%                         b(e) = im_s(y,x,c) - im_s(y,x-1,c) + im_background(y,x-1,c);

                        s = im_s(y,x,c) - im_s(y,x-1,c);
                        t = im_background(y,x,c) - im_background(y,x-1,c);
                        d = alfa*s + (1-alfa)*t;
                        
                        b(e) = d + im_background(y,x-1,c);

                    end
                    e=e+1;
                end
                X=X+1;
            end
            Y=Y+1;
        end

        v{c} = A\b;
    end
    
    for c = 1:3 
       im_out(:,:,c)  =  reshape(v{c}, [imh imw]);
    end

    im_final = im_background;
    im_final( j_beg:j_end, k_beg:k_end, : ) = im_out;

end