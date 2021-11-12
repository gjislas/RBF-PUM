%%
clc
clear

N = 100000;

TestD = rand(N);
TestS = single(TestD);

wBinary(TestD,'100kx100k_D','double');
wBinary(TestS,'100kx100k_S','single');