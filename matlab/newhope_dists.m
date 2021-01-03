eta = 8;

n = 1024;
tests = 10000; % gets squared!

cbd1 = binornd(2*eta, 0.5,  [tests, n]) - eta;
cbd2 = binornd(2*eta, 0.5, [n, tests]) - eta;

res = cbd1*cbd2;
clear cbd1 cbd2

figure; histogram(res(:));
