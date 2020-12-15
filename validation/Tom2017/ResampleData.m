% resample data

data = readtable('.//Tom2017.csv');

wstart  = max(min(data.wstar_mu),min(data.wstar_nu))
wend    = min(max(data.wstar_mu),max(data.wstar_nu)) % maybe min

x1 = data.wstar_mu(1:end-1)
v1 = data.mu55star(1:end-1)
xq = linspace(wstart,wend,100)
vq1 = interp1(x1,v1,xq);



x2 = data.wstar_nu
x2 = x2(~isnan(x2))
v2 = data.nu55star
v2 = v2(~isnan(x2))
xq = linspace(wstart,wend,100)
vq2 = interp1(x2,v2,xq);

wstar    = xq.'
mu55star = vq1.'
nu55star = vq2.'

T = table(wstar,mu55star,nu55star)

hold on
plot(xq,vq1)
plot(xq,vq2)


writetable(T,'Tom2017resample.csv')