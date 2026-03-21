function s=UU(x, sigma)

if sigma>=0.99 && sigma<= 1.01
    s=log(x);
else
    s=(x.^(1-sigma) -1)/(1-sigma);
end

end