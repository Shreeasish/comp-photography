function kernel = compute_kernel(sigma,size)
[x,y] = meshgrid(-size/2:size/2,-size/2:size/2);
constant = 1/(2*pi*sigma*sigma);
kernel = constant*exp( -(y.^2 + x.^2 )/(2 * sigma * sigma));

kernel=kernel/sum(kernel(:));
end
% kernel = kernel/sum(kernel(:));
% imagesc(kernel);
