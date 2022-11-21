%%NULLCLINES
clear all;
close all hidden;
clc;

syms M H R a1 a2 a3 a4 a5 a6 a7;
f1 = 1+a1*M*(1-M)-a2*M*H;
f2 = a3*H*R - a4*H;
f3 = a5*R*(1-R)-a6*H*R-a7*R;
cond = [a1>0 a2>0 a3>0 a4>0 a5>0 a6>0 a7>0 M>=0 H>=0 R>=0];
assume(cond);

M1 = solve(f1==0,M);
H1 = solve(f1==0,H);
H2 = solve(f2==0,H);
H3 = solve(f3==0,H);
R1 = solve(f2==0,R);
R2 = solve(f3==0,R);

NM = [M1];
NH = [H1; H2; H3];
NR = [R1; R2];
disp('NM');
disp(NM);
disp('NH');
disp(NH);
disp('NR');
disp(NR);

%%
%%FIXEDPOINTS
clear all;
close all hidden;
clc;
tic
syms M H R a1 a2 a3 a4 a5 a6 a7;
f1 = 1+a1*M*(1-M)-a2*M*H;
f2 = a3*H*R - a4*H;
f3 = a5*R*(1-R)-a6*H*R-a7*R;
%cond = [a1>=0 a2>=0 a3>=0 a4>=0 a5>=0 a6>=0 a7>=0 M>=0 H>=0 R>=0]; % 32
%fixed points!
%cond = [a1==0 a2==0 a3==0 a4==0 a5==0 a6==0 a7==0 M==0 H==0 R==0]; % no
%fixed points!
cond = [a1>0 a2>0 a3>0 a4>0 a5>0 a6>0 a7>0 M>=0 H>=0 R>=0]; % 6 fixed points!
assume(cond);
S = solve([f1==0 f2==0 f3==0],[M H R],'Real',true,'Returnconditions',true);
toc
Fixedpoints = simplify([S.M S.H S.R]);
Parameters = simplify(S.parameters);
Conditions = simplify(S.conditions);
toc
disp('Fixedpoints');
pretty((Fixedpoints));
simplify(Fixedpoints)
disp('Conditions');
pretty((Conditions));
disp('Parameters');
pretty((Parameters));
toc

%%
%%NUMERIC FIXED POINTS
clear all;
close all hidden;
clc;

condition = 1;
switch condition
    case 1
        a = [4.5 1 0.5 0.4 0.7 0.2 0.4];
    case 2
        a = [3.5 1 5 0.4 0.7 0.1 0.1];
    case 3
        a = [3 1 4.8 0.4 3.7 1.9 0.1];
    otherwise
        disp('Invalid condition, choose 1, 2 or 3.');
        return;
end

syms M H R;

f1 = 1+a(1)*M*(1-M)-a(2)*M*H;
f2 = a(3)*H*R - a(4)*H;
f3 = a(5)*R*(1-R)-a(6)*H*R-a(7)*R;

assume([M>=0 H>=0 R>=0]);
S = vpasolve([f1==0 f2==0 f3==0],[M H R]);
S.M;
S.H;
S.R;

Fixedpoints = [double(S.M) double(S.H) double(S.R)];

stability = [];
for s = 1:size(Fixedpoints,1)
M0 = Fixedpoints(s,1);
H0 = Fixedpoints(s,2);
R0 = Fixedpoints(s,3);

A = [a(1)-2*a(1)*M0-a(2)*H0 -a(2)*M0     0;
     0                      a(3)*R0-a(4) a(3)*H0;
     0                      -a(6)*R0     a(5)-2*a(5)*R0-a(6)*H0-a(7)];
 
E = eig(A);
disp(E);
if(max(real(E)) < 0)
    stability = [stability; 1];
else
    stability = [stability; 0];
end
end

Fixedpoints = [Fixedpoints stability];

disp('Fixedpoints');
disp('|    m    |    h    |    r    |stability|');
disp(Fixedpoints);
