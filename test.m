clear all
close all
clc

t = 1:500;
c = 2.8;

alpha = c./(1+exp(-1e5*(t-100)));

plot(alpha)