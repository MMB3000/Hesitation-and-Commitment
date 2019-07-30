function [ gn, gpn, dgn, dgpn, pggn ] = new_new_belief(alfa,nmax)
alpha = alfa;
if nargin<2
    nmax = 28;
end
gn = {};
gpn = {};
dgn = {};
dgpn = {};
pggn = {};


for n = 1:nmax
    r   = linspace(0,n,n+1);
    g_numerator = bsxfun(@minus,betainc(1-alpha,r+1,n-r+1),betainc(0.5,r+1,n-r+1));
    g_denominator = bsxfun(@minus,betainc(1-alpha,r+1,n-r+1),betainc(alpha,r+1,n-r+1));
    g = bsxfun(@rdivide,g_numerator,g_denominator);    
    
    dg  = g(1:end-1)+diff(g)/2;
    dg  = diff(dg);
    
    r_n = linspace(0,n+1,n+2);
    g_n_numerator = bsxfun(@minus,betainc(1-alpha,r_n+1,n-r_n+2),betainc(0.5,r_n+1,n-r_n+2));
    g_n_denominator = bsxfun(@minus,betainc(1-alpha,r_n+1,n-r_n+2),betainc(alpha,r_n+1,n-r_n+2));
    g_n = bsxfun(@rdivide,g_n_numerator,g_n_denominator);
    
    dg_n= g_n(1:end-1)+diff(g_n)/2;
    dg_n= diff(dg_n);
    
    gn{end+1} = g';
    gpn{end+1} = g_n';
    dgn{end+1} = dg';
    dgpn{end+1} = dg_n';
    
    pgg         = compute_pgg(r,n,alpha,r_n);
    pggn{end+1} = pgg;
        
    
end

end

function pgg = compute_pgg(r,n,alfa,r_n)

for i = 1:size(r,2)
    for j = 1:size(r_n,2)
        if i<=j && j<=i+1
            pgg(i,j) = beta(r_n(1,j)+1,n+1-r_n(1,j)+1)/beta(r(1,i)+1,n-r(1,i)+1)*(betainc(1-alfa,r_n(1,j)+1,n+1-r_n(1,j)+1)-betainc(alfa,r_n(1,j)+1,n+1-r_n(1,j)+1))/(betainc(1-alfa,r(1,i)+1,n-r(1,i)+1)-betainc(alfa,r(1,i)+1,n-r(1,i)+1));
        else
            pgg(i,j) = 0;
        end
    end
end

end