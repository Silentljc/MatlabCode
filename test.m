clear;clc;
a=ones(7,1);
b=ones(8,1);
c=ones(9,1)*4;
d=ones(8,1);
e=ones(7,1);
f=[3;4;ones(5,1)*5;4;3]+ones(9,1)*3;

x = chase5(a,b,c,d,e,f)

err = ones(9,1)-x
 
