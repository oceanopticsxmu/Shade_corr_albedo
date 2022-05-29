
function [shade]=get_shderror(a,bb,solz)

% calculate shading error from IOPs and solz
s0(1:6)=[0.0137,0.0759,0.1328,0.2298,5,0.6874];
t1=(s0(1).*log(90-solz)+s0(2)).*exp(-0.1.*a);
t2=(s0(3).*log(90-solz)+s0(4)).*(a+bb).*exp(-s0(5)*bb);
t3=s0(6).*a.*bb;

shade=1-exp(-(t1+t2+t3));

end