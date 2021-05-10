function dpbdt = recdefODEs(t,pb,k_r,k_d,m,n)
% Code for equations (2.3) & (2.4)

% pb(1) = p;
% pb(2) = b;

dpbdt = [k_r.*(pb(1)^m).*pb(2) - k_d.*pb(1).*(pb(2)^n);
            k_d.*pb(1).*(pb(2)^n) - k_r.*(pb(1)^m).*pb(2)];
end