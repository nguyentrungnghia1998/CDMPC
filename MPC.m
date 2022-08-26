%% Simulation DMPC
clc;
clear;
close all;
%% Time and step
Step = 0.2;
T_end =  50;
t = 0:Step:T_end;
%% Variable
y1d = cell(1,length(t));
y2d = cell(1,length(t));
x1 = cell(1,length(t));
x2 = cell(1,length(t));
x_1pred = cell(1,length(t));
x_2pred = cell(1,length(t));
u1 = cell(1,length(t));
u2 = cell(1,length(t));
u_1pred = cell(1,length(t));
u_2pred = cell(1,length(t));
y1 = cell(1,length(t));
y2 = cell(1,length(t));
%% Parameters
alpha = 10;
gamma = 1;
P = 100;
M = P;
A11 = [2.74 -1.27 0.79 0;
       2 0 0 0;
       0 0.5 0 0;
       0 0 0 0.37];
A22 = [1.68 -0.82 0 0;
       1 0 0 0;
       0 0 1.57 -0.67;
       0 0 1 0];
A = [A11 zeros(4);
     zeros(4) A22];
B1 = [0.25;0;0;0.5;0;0;0;0];
B2 = [0;0;0;0;0.25;0;0.5;0];
B11 = [0.25;0;0;0.5];
B22 = [0.25;0;0.5;0];
B = [B1 B2];
C = [-0.1 0.03 0.12 0 0.07*alpha 0.07*alpha 0 0;
     0 0 0 0.25*alpha 0 0 0.29 -0.2];
C11 = [-0.1 0.03 0.12 0];
C12 = alpha*[0.07 0.07 0 0];
C21 = alpha*[0 0 0 0.25];
C22 = [0 0 0.29 -0.2];
Q = [1 0;0 1];
Q_ = eye(2*M);
R1_ = gamma*eye(M);
R2_ = gamma*eye(M);
GAMMA1 = [1 zeros(1,M-1)];
GAMMA2 = GAMMA1;
GAMMA_1 = tril(ones(M));
GAMMA_2 = tril(ones(M));
C_cell = repmat({C},1,P);
Ca_ = blkdiag(C_cell{:});
S = cell(P);
for i1 = 1:P
    for j1=1:P
        if i1-j1<0
            S{i1,j1}=zeros(8);
        else
            S{i1,j1}=A^(i1-j1);
        end
    end
end
S_ = cell2mat(S);
S = blkdiag(S_,S_);
B1_cell = repmat({B1},1,M-1);
B2_cell = repmat({B2},1,M-1);
B1_ = [zeros(8*(M-1),1) blkdiag(B1_cell{:});
       zeros(8,(M-1)*1) B1];
B2_ =  [zeros(8*(M-1),1) blkdiag(B2_cell{:});
       zeros(8,(M-1)*1) B2];
Aa_ = [A;zeros(8*(M-1),8)];
L1 = [eye(4);zeros(4)];
L2 = [zeros(4);eye(4)];
L1_dot = [zeros(4) zeros(4);zeros(4) eye(4)];
L2_dot = [eye(4) zeros(4);zeros(4) zeros(4)];
T_nga1 = [0 0;0 1];
T_nga2 = [1 0;0 0];
BT_1_cell = repmat({B*T_nga1},1,M-1);
BT_2_cell = repmat({B*T_nga2},1,M-1);
B1_nga = [zeros(8*(M-1),2) blkdiag(BT_1_cell{:});
       zeros(8,(M-1)*2) B*T_nga1];
B2_nga = [zeros(8*(M-1),2) blkdiag(BT_2_cell{:});
       zeros(8,(M-1)*2) B*T_nga2];
N1 = Ca_*S_*B1_*GAMMA_1;
N2 = Ca_*S_*B2_*GAMMA_2;
H1 = N1'*Q_*N1+R1_;
H2 = N2'*Q_*N2+R2_;
K1_ = pinv(H1)*N1'*Q_;
K2_ = pinv(H2)*N2'*Q_;
K1 = GAMMA1*K1_;
K2 = GAMMA2*K2_;
GAMMA_1_dot = ones(M,1);
GAMMA_2_dot = ones(M,1);
%% Initial value
x1{1} = [0;0;0;0];
x2{1} = [0;0;0;0];
u1{1} = 0;
u2{1} = 0;
u1_pred_old = zeros(M,1);
u2_pred_old = zeros(M,1);
x1_pred_old = zeros(8,1);
x2_pred_old = zeros(8,1);
U_pred_old = reshape([u1_pred_old u2_pred_old]',[2*M 1]);
%% Simulation
for i = 1:length(t)
    if t(i)<10
        y1d{i} = 0;
    elseif t(i)<20
        y1d{i} = 0.5;
    elseif t(i)<25
        y1d{i} = 0.75;
    else
        y1d{i} = 0;
    end
    if t(i)<5
        y2d{i} = 0;
    elseif t(i)<20
        y2d{i} = 0.5;
    elseif t(i)<25
        y2d{i} = 0.25;
    elseif t(i)<35
        y2d{i} = 0;
    else 
        y2d{i} = 0.5;
    end
end

y1d = cell2mat(y1d);
y2d = cell2mat(y2d);
y1d_copy = y1d;
y2d_copy = y2d;
y1d(end+1:end+P) = y1d(end);
y2d(end+1:end+P) = y2d(end);
yd = [y1d;y2d];
for i = 1:length(t)
    y1{i} = C11*x1{i}+C12*x2{i};
    y2{i} = C21*x1{i}+C22*x2{i};
    if i~=1
        Yd = yd(:,i+1:i+P);
        Yd = reshape(Yd,2*P,1);
        Z1m = Ca_*S_*(B1_*GAMMA_1_dot*u1{i-1}+Aa_*L1*x1{i}+Aa_*L1_dot*x1_pred_old+B1_nga*U_pred_old);
        Z2m = Ca_*S_*(B2_*GAMMA_2_dot*u2{i-1}+Aa_*L2*x2{i}+Aa_*L2_dot*x2_pred_old+B2_nga*U_pred_old);
        u1{i} = u1{i-1} + K1*(Yd-Z1m);
        u2{i} = u2{i-1} + K2*(Yd-Z2m);
        u1_pred_new = GAMMA_1_dot*u1{i-1} + K1_*(Yd-Z1m);
        u2_pred_new = GAMMA_2_dot*u2{i-1} + K2_*(Yd-Z2m);
        x1_pred_new = A*L1*x1{i} + A*L1_dot*x1_pred_old + B1*u1_pred_new(2)+B2*u2_pred_old(3);
        x2_pred_new = A*L2*x2{i} + A*L2_dot*x2_pred_old + B2*u2_pred_new(2)+ B1*u1_pred_old(3);
        U_pred_new = reshape([u1_pred_new u2_pred_new]',[2*M 1]);
        u1_pred_old = u1_pred_new;
        u2_pred_old = u2_pred_new;
        x1_pred_old = x1_pred_new;
        x2_pred_old = x2_pred_new;
        U_pred_old = U_pred_new;
    end
    if i == length(t)
        break
    end
    %% Update
    x1{i+1} = A11*x1{i} + B11*u1{i};
    x2{i+1} = A22*x1{i} + B22*u2{i};
end

y1 = cell2mat(y1);
y2 = cell2mat(y2);
figure(1);
plot(t,y1d_copy,t,y1);
figure(2);
plot(t,y2d_copy,t,y2);