clear;
clc;
close all;
Ps = 2:2:30;
log_gamma = -3:0.2:3;
Ni = length(log_gamma);
Nj = length(Ps);
max_eig_value = zeros(Nj,Ni);
alpha = 0.1;
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
B = [B1 B2];
L1 = [eye(4);zeros(4)];
L2 = [zeros(4);eye(4)];
L = blkdiag(L1,L2);
L1_dot = [zeros(4) zeros(4);zeros(4) eye(4)];
L2_dot = [eye(4) zeros(4);zeros(4) zeros(4)];
L_dot = [L1_dot;L2_dot];
Q = [1 0;0 1];
C = [-0.1 0.03 0.12 0 0.07*alpha 0.07*alpha 0 0;
     0 0 0 0.25*alpha 0 0 0.29 -0.2];
for i = 1:Ni
    for j = 1:Nj
        gamma = 10^log_gamma(i);
        P = Ps(j);
        M = P;
        AN = cell(4,4);
        AN{1,1} = A;
        GAMMA1 = [1 zeros(1,M-1)];
        GAMMA2 = [1 zeros(1,M-1)];
        GAMMA = blkdiag(GAMMA1,GAMMA2);
        AN{1,3} = B*GAMMA;
        L1_cell = repmat({L1'},1,P);
        L2_cell = repmat({L2'},1,P);
        L1_ = blkdiag(L1_cell{:});
        L2_ = blkdiag(L2_cell{:});
        L_ = blkdiag(L1_,L2_);
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
        Aa_ = [A;zeros(8*(M-1),8)];
        A_ = blkdiag(Aa_,Aa_);
        AN{2,1} = L_*S*A_*L;
        L_nga = L_dot * [eye(8) zeros(8,8*(P-1))];
        OMEGA = cell(P,1);
        OMEGA11 = [eye(4) zeros(4,4*(P-1))];
        for i1 = 1:P
            tg = circshift(OMEGA11,4*(i1-1),2);
            OMEGA{i1} =blkdiag(tg,tg);
        end
        OMEGA = cell2mat(OMEGA);
        AN{2,2} = L_*S*A_*L_nga*OMEGA;
        B1_cell = repmat({B1},1,M-1);
        B2_cell = repmat({B2},1,M-1);
        B1_ = [zeros(8*(M-1),1) blkdiag(B1_cell{:});
               zeros(8,(M-1)*1) B1];
        B2_ =  [zeros(8*(M-1),1) blkdiag(B2_cell{:});
               zeros(8,(M-1)*1) B2];
        B_ = blkdiag(B1_,B2_);
        AN{2,3} = L_*S*B_;
        T_nga1 = [0 0;0 1];
        T_nga2 = [1 0;0 0];
        BT_1_cell = repmat({B*T_nga1},1,M-1);
        BT_2_cell = repmat({B*T_nga2},1,M-1);
        B1_nga = [zeros(8*(M-1),2) blkdiag(BT_1_cell{:});
               zeros(8,(M-1)*2) B*T_nga1];
        B2_nga = [zeros(8*(M-1),2) blkdiag(BT_2_cell{:});
               zeros(8,(M-1)*2) B*T_nga2];
        B_nga = [B1_nga;B2_nga];
        PI = cell(M,1);
        PI11 = [eye(1) zeros(1,1*(M-1))];
        for i1 = 1:M
            tg = circshift(PI11,1*(i1-1),2);
            PI{i1} =blkdiag(tg,tg);
        end
        PI = cell2mat(PI);
        AN{2,4} = L_*S*B_nga*PI;
        GAMMA_1 = tril(ones(M));
        GAMMA_2 = tril(ones(M));
        Q_ = eye(2*M);
        C_cell = repmat({C},1,P);
        Ca_ = blkdiag(C_cell{:});
        C_ = blkdiag(Ca_,Ca_);
        N1 = Ca_*S_*B1_*GAMMA_1;
        N2 = Ca_*S_*B2_*GAMMA_2;
        R1_ = gamma*eye(M);
        R2_ = gamma*eye(M);
        H1 = N1'*Q_*N1 + R1_;
        H2 = N2'*Q_*N2 + R2_;
        K1_ = pinv(H1)*N1'*Q_;
        K2_ = pinv(H2)*N2'*Q_;
        PSI = blkdiag(GAMMA_1*K1_,GAMMA_2*K2_);
        Os = -PSI*C_*S*A_*L;
        THETA = -PSI*C_*S*A_*L_nga*OMEGA;
        GAMMA_1_dot = ones(M,1);
        GAMMA_2_dot = ones(M,1);
        GAMMA_dot = blkdiag(GAMMA_1_dot,GAMMA_2_dot);
        Yhoa = GAMMA_dot*GAMMA - PSI*C_*S*(B_*GAMMA_dot*GAMMA+B_nga*PI);
        AN{3,1} = Os*A+THETA*L_*S*A_*L;
        AN{3,2} = THETA*L_*S*A_*L_nga*OMEGA;
        AN{3,3} = Os*B*GAMMA + THETA*L_*S*B_ + Yhoa;
        AN{3,4} = THETA*L_*S*B_nga*PI;
        AN{4,3} = eye(2*M);
        AN{1,2} = zeros(8,8*P);
        AN{1,4} = zeros(8,M*2);
        AN{4,1} = zeros(2*M,8);
        AN{4,2} = zeros(2*M,8*P);
        AN{4,4} = zeros(2*M,2*M);
        AN_copy = AN;
        AN = cell2mat(AN);
        max_eig_value(j,i) = max(abs(eig(AN)));
    end
end