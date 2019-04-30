fsamp = 6.28;
fcuts = [0.78 1.31 1.83 2.36];
mags = [1 0 1];
devs = [0.001 0.001 0.001];

[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
n = n + rem(n,2);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');

[H,f] = freqz(hh,1,1024,fsamp);
plot(f,abs(H))
grid