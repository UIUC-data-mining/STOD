function x = ktensor3(vec,mat)
% Chi Wang's own implementation of tensor(ktensor(vec,mat,mat,mat))
n = size(mat,2);
m = size(mat,1);
x = zeros(m*m*m,1);
for i=1:n
    mi = mat(:,i);
    p = mi * mi';
    p = p(:) * mi';
    x = x + vec(i)*p(:);
end
x = reshape(x,[m m m]);
