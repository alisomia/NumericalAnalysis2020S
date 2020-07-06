clf;
f = @(t) exp(-t.*t/10).*(sin(2*t)+2*cos(4*t)+0.4*sin(t).*sin(50*t));
X = 0:1/256:255/256;
X = X*pi;
plot(X,f(X));
hold on
out = textread('./out50');
plot(X,out)

