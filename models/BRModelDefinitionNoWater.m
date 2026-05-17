% SCRIPT WITH SOME MODEL DESCRIPTIONS TO TEST SIMULATION ALGORITHMS
%If we don't consider the water we have the following, with the convention
%       1   2    3    4  5   6    7   8   9    10   11   12  13    14
%       IO3 H2O2 HO2  H  Mn  MnOH O2  MA IMA  IO2  HIO2  I   HOI   I2
%R1:    IO3+I+2*H -k1-> HIO2+HOI
%R2:    HIO2+I+H -k2->H 2HOI
%R3:    HOI+I+H -k3-> 2HOI
%R4:    IO3+HIO2 + H -k4-> 2IO2
%R5:    2*HIO2 -k5-> IO3+HOI+H
%R6:    O2+Mn - c*k6 ->HIO2+MnOh
%R7:    MnOH+H2O2 - k7 -> Mn+HO2
%R8:    2*HO2 +H2O2 - k8 -> H2O2+O2     
%R9:    I2+MA - k9 -> IMA + I + H
%R10:   HOI+H2O2 - k10 -> I + O2 +H
%R-3:   I2 - c*k-3 - > HOI + I + H
%R-4:   2IO2 - c* k-4 -> IO3+HIO2
avogadro = 6.02214076e23;

vMinusBR =        [1 0 0 2 0 0 0 0 0 0 0 1 0 0; %1
                   0 0 0 1 0 0 0 0 0 0 1 1 0 0; %2
                   0 0 0 1 0 0 0 0 0 0 0 1 1 0; %3
                   1 0 0 1 0 0 0 0 0 0 1 0 0 0; %4
                   0 0 0 0 0 0 0 0 0 0 2 0 0 0; %5
                   0 0 0 0 1 0 0 0 0 1 0 0 0 0; %6
                   0 1 0 0 0 1 0 0 0 0 0 0 0 0; %7
                   0 0 2 0 0 0 0 0 0 0 0 0 0 0; %8
                   0 0 0 0 0 0 0 1 0 0 0 0 0 1; %9
                   0 1 0 0 0 0 0 0 0 0 0 0 1 0;%10
                   0 0 0 0 0 0 0 0 0 0 0 0 0 1; %-3
                   0 0 0 0 0 0 0 0 0 2 0 0 0 0]; %-4
                   
vPlusBR   =        [0 0 0 0 0 0 0 0 0 0 1 0 1 0; %1
                    0 0 0 0 0 0 0 0 0 0 0 0 2 0; %2
                    0 0 0 0 0 0 0 0 0 0 0 0 0 1; %3
                    0 0 0 0 0 0 0 0 0 2 0 0 0 0; %4
                    1 0 0 1 0 0 0 0 0 0 0 0 1 0; %5
                    0 0 0 0 0 1 0 0 0 0 1 0 0 0; %6
                    0 0 1 0 1 0 0 0 0 0 0 0 0 0; %7
                    0 1 0 0 0 0 1 0 0 0 0 0 0 0; %8
                    0 0 0 1 0 0 0 0 1 0 0 1 0 0; %9 
                    0 0 0 1 0 0 1 0 0 0 0 1 0 0; %10
                    0 0 0 1 0 0 0 0 0 0 0 1 1 0; %-3
                    1 0 0 1 0 0 0 0 0 0 1 0 0 0]; %-4
vBR = vPlusBR - vMinusBR;
c=55.5;
a = 0.0225;
b = 0.33;
h = 0.056;
q = 0.0015;
c = 55.5; 
v0 = 0.004;
kBR = [5*1e3 2*1e10 3.1*1e12 7.35*1e3 1e5 c*180.2 3.2*1e4 7.5*1e5 12.5 37 c*0.013 c*8*1e5];
[cBR,orderBR]=generateStochasticrates(kBR,vMinusBR,1e8);

initialMolarConcentrationsBR = [a b 1e-3 h 0.056 v0 0 q 0 1e-5 0 4e-4 0 0]; 


%These causes the oscillations of HIO, HIO2
%initialStateBR = 2*ones(1, 14)*1000;
%initialStateBR(2)=1e4;
%initialStateBR(11)=0;
%initialState(14)=0;
%initialState(6)=0;
%MNOH=0; HIO2=0 , IO2=0; 
initialStateBR=initialMolarConcentrationsBR*1e8;


%10 10 10 11 10 10 2245
%"I","HOI","HIO2","I2","IO2","Mn","HO2"