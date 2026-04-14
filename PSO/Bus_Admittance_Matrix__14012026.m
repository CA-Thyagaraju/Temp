clc;
clear;
y = input('Enter the admittance matrix y: ');
Ybus = formYbus(y);
yaa = ybus(y);
fprintf('Bus Admittance Matrix Ybus =\n');
disp(Ybus)
disp(yaa)