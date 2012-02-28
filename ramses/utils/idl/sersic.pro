function sersic,r,n=n,re=re,se=se

k=2.*n-0.327
return,se*exp(-k*((r/re)^(1./n)-1))

end
