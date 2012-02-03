function A = gaussTransform(in, mu, Sigma)
mu_row = mu(:)';
R = chol(Sigma);
A = repmat(mu,length(in),1)+in*R;
end
