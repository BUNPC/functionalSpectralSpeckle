function PathLengths = pathlengths(lambda)

Hb0 = 60*10^-6;
HbR = 40*10^-6;
% lambda = [470, 530, 590, 625];
g = 0.9;
c = 3*10^10;
e = GetExtinctions(lambda);
mua = e(:,1).*Hb0+e(:,2).*HbR;
mus = 150*(lambda/560).^(-2);
z0 = 1./((1-g)*mus');
gamma = sqrt(c./(3*(mua+(1-g)*mus')));
PathLengths = (c*z0./(2*gamma.*(mua*c).^0.5)).*(1+(3/c)*mua.*(gamma).^2);