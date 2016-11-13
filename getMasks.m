function [masks] = getMasks(mu,r1,r2,alpha)

mu = mu';

w = complex(cos(2*pi/3),sin(2*pi/3));
M = length(mu);

for h = 1:M
    muAlpha(1,h) = exp(2*pi*sqrt(-1)*alpha*(h-1)/M) * mu(r2,h);
end

masks = [(mu(r1,:)+w^0*muAlpha);...
    (mu(r1,:)+w^1*muAlpha);...
    (mu(r1,:)+w^2*muAlpha)];

end