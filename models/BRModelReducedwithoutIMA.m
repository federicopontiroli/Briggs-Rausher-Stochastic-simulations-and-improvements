% SCRIPT WITH SOME MODEL DESCRIPTIONS TO TEST SIMULATION ALGORITHMS
%If we don't consider the water we have the following, with the convention
%       1   2    3    4  5    6    7       8 
%       I HOI  HIO2  I2  IO2  Mn   HO2     MnOH
vMinusBRr =      [ 1 0 0 0 0 0 0 0   %1
                   1 0 1 0 0 0 0 0 ; %2
                   1 1 0 0 0 0 0 0 ; %3
                   0 0 1 0 0 0 0 0 ; %4
                   0 0 2 0 0 0 0 0; %5
                   0 0 0 0 1 1 0 0 ; %6
                   0 0 0 0 0 0 0 1 ; %7
                   0 0 0 0 0 0 2 0 ; %8
                   0 0 0 1 0 0 0 0 ; %9
                   0 1 0 0 0 0 0 0 ;%10
                   0 0 0 1 0 0 0 0 ; %-3
                   0 0 0 0 2 0 0 0 ]; %-4
                   
vPlusBRr   =      [ 0 1 1 0 0 0 0 0 ; %1
                    0 2 0 0 0 0 0 0 ; %2
                    0 0 0 1 0 0 0 0 ; %3
                    0 0 0 0 2 0 0 0 ; %4
                    0 1 0 0 0 0 0 0 ; %5
                    0 0 1 0 0 0 0 1 ; %6
                    0 0 0 0 0 1 1 0 ; %7 
                    0 0 0 0 0 0 0 0 ; %8
                    1 0 0 0 0 0 0 0; %9
                    1 0 0 0 0 0 0 0; %10
                    1 1 0 0 0 0 0 0; %-3
                    0 0 1 0 0 0 0 0 ]; %-4

vBRr = vPlusBRr - vMinusBRr;
a = 0.0225;
b = 0.33;
h = 0.056;
q = 0.0015;
c = 55.5; 
v0 = 0.004;

kBR = [a*(h^2)*5*1e3 h*2*1e10 h*3.1*1e12 a*h*7.35*1e3 1e5 c*180.2 b*3.2*1e4 7.5*1e5 q*12.5 b*37 c*0.013 c*8*1e5];
[cBRr,orderBR]=generateStochasticrates(kBR,vMinusBRr,1e8);
%y0=[1e-5 0 0 1e-5 1e-5 v0 10^(-3) 1e-5 1e-5 0];
y0=[0.0004 0 0 0 1e-5 v0 10^(-3) 0];

initialStateBRr = round(1e8 * y0);
