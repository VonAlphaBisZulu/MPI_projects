load('duality_model.mat');
cnap = model;
cnap.reacMin(ismember(cellstr(cnap.reacID),'r_6')) = 0;
cnap.reacMax(ismember(cellstr(cnap.reacID),'r_6')) = 0;
% cnap.reacMin(ismember(cellstr(cnap.reacID),'r_7')) = 0;
% cnap.reacMax(ismember(cellstr(cnap.reacID),'r_7')) = 0;
% cnap.reacMin(ismember(cellstr(cnap.reacID),'r_4')) = 0; % needed for essential production
% cnap.reacMax(ismember(cellstr(cnap.reacID),'r_4')) = 0; % needed for essential production
% translating
A  = cnap.stoichMat;
b  = zeros(cnap.nums,1);
lb = cnap.reacMin;
ub = cnap.reacMax;
c  = cnap.objFunc;
n  = size(cnap.stoichMat,2);

%% 0. Test original model with linprog
linprog_options = optimoptions('linprog','Display','off');
intlinprog_options = optimoptions('intlinprog','Display','off');

% opt_p0 is the optimal value in the primal problem
[~,opt_ref,~] = linprog(c,[],[],A,b,lb,ub,[],linprog_options);



%% 1. Split all reactions, optimize with linprog and verify if identical solution is found
A_p1  = [A  , -A ];
ub1 = max(zeros(2*n,1),[ub;-lb]);
lb1 = max(zeros(2*n,1),[lb;-ub]);
c_p1 = [c;-c];

[~,opt_p1,~] = linprog(c_p1,[],[],A_p1,b,lb1,ub1,[],linprog_options);

if isempty(opt_p1) || abs(opt_ref-opt_p1)>1e-7
    error('LP didn''t return correct solution');
end

%% 2. Reshape this problem into A*x<=b, x>=0 and min(c'*x)
A_p2 = [A_p1 ;-A_p1;eye(2*n);-eye(2*n)]; % to mime A*x = 0 and upper bound
b_p2 = [ b ; b ;  ub1   ;  -lb1   ];
c_p2 = c_p1;

[~,opt_p2,~] = linprog(c_p2,A_p2,b_p2,[],[],zeros(length(c_p2),1),inf(length(c_p2),1),[],linprog_options);

if isempty(opt_p2) || abs(opt_ref-opt_p2)>1e-7
    error('LP didn''t return correct solution');
end

%% 3. Solve dual problem
%
% Primal: -min(-c'x),  A *x<=b, x>=0
% Dual:    min(b'y),   A'*y>=c, y>=0
% with:    min(c'x) = max(b'y);
%
% => DUAL: min(b'y), -A'*y<=-c, y>=0 
%
%
A_d1 = -A_p2';
b_d1 = c_p2;
c_d1 = b_p2;
% 
[~,q_d1,~] = linprog(c_d1,A_d1,b_d1,[],[],zeros(length(c_d1),1),inf(length(c_d1),1),[],linprog_options);

if isempty(q_d1) || abs(opt_ref+q_d1)>1e-7
    error('LP didn''t return correct solution');
end

%% 4. Connect to primal problem - bilevel optimization
A_ineq_bi1 = [zeros(size(A_d1,1),size(A,2)) , A_d1];
b_ineq_bi1 = b_d1;
A_eq_bi1 = [A , zeros(size(A,1),size(A_d1,2)) ];
b_eq_bi1 = b;
% adding a constraint to connect both biomass rates
A_eq_bi1 = [A_eq_bi1 ; c' c_d1'];
b_eq_bi1 = [b_eq_bi1; 0];

lb_bi1 = [lb; zeros(length(c_d1),1)];
ub_bi1 = [ub; inf(length(c_d1),1)];

c_prod_bi1 = zeros(1,size(A_ineq_bi1,2)); % objective function is product synthesis
c_prod_bi1(2) = -1;
[~,q_d2,~] = linprog(c_prod_bi1, A_ineq_bi1, b_ineq_bi1, A_eq_bi1, b_eq_bi1, lb_bi1, ub_bi1,[],linprog_options);

%% 5. Primal problem with knockable reaction 6 / MILP

idx_r6 = find(ismember(cellstr(cnap.reacID),'r_6'));
cnap.reacMax(idx_r6) = 100;

A  = cnap.stoichMat;
b  = zeros(cnap.nums,1);
lb = cnap.reacMin;
ub = cnap.reacMax;
c  = cnap.objFunc;
n  = size(cnap.stoichMat,2);
% 5.1 Split all reactions
A_pk  = [A  , -A ];
ub1 = max(zeros(2*n,1),[ub;-lb]);
lb1 = max(zeros(2*n,1),[lb;-ub]);
c_p1 = [c;-c];
% 5.2 Reshape into A*x<=b, x>=0 and min(c'*x)
A_pk2 = [A_pk ;-A_pk;eye(2*n);-eye(2*n)]; % to mime A*x = 0 and upper bound
b_pk2 = [ b ; b ;  ub1   ;  -lb1   ];
c_pk2 = c_p1;
% 5.3 Adapt lb/ub rows for reaction 6 to account for reaction-switch and add variable
A_pk3 = [A_pk2 , zeros(size(A_pk2,1),1)];
A_pk3(2*size(A_pk,1)+            idx_r6,end) =  ub1(idx_r6);
A_pk3(2*size(A_pk,1)+size(A,2)  +idx_r6,end) =  ub1(size(A,2)+idx_r6);
A_pk3(2*size(A_pk,1)+2*size(A,2)+idx_r6,end) = -lb1(idx_r6);
A_pk3(2*size(A_pk,1)+3*size(A,2)+idx_r6,end) = -lb1(size(A,2)+idx_r6);
b_pk3 = b_pk2;
c_pk3 = [c_pk2 ; 0];
% 5.4 Add two rows to constrain the bool variable z <= 1 and z >= 1 or 0. (-z <= -1)
A_pk4 = [A_pk3; zeros(2,size(A_pk3,2))];
A_pk4((end-1:end),end) = [1 -1];
c_pk4 = c_pk3;
b_pk4 = [b_pk3; 1; -1];
lb_pk4 = zeros(size(A_pk4,2),1); % reaction open or closed possible
ub_pk4 = inf(size(A_pk4,2),1); 
ctype_pk4 = [repmat('C',1,2*size(A,2)),'B'];

% with reaction 6 off (z >= 1) should restrict growth to 3.5
[~,opt_pk1,~] = intlinprog(c_pk4, find(ctype_pk4=='B'), A_pk4, b_pk4, [], [], lb_pk4, ub_pk4,[],intlinprog_options);

% with reaction 6 on (z >= 0) should permit a higher maximum growth (7)
b_pk4(end) = 0;
[~,opt_pk2,~] = intlinprog(c_pk4, find(ctype_pk4=='B'), A_pk4, b_pk4, [], [], lb_pk4, ub_pk4,[],intlinprog_options);

if -opt_pk1 >= -opt_pk2
    error('LP didn''t return correct solution');
end

%% Hier habe ich die Dualitäts-Funktion aus RobustKnock verwendet.
%  Die Notation ist anders, da z = 1 heißt, dass die Reaktion NICHT ausgeknockt ist
%  Im Grunde egal wierum es funktioniert, wenn wir eine Richtung hinkriegen, dann auch die Andere.
% Ansonsten:
% wir -> RobustKnock
%  A -> [A_w Ay_w]
%  x -> v
%  y -> w
%  z -> y (invers)
%  p -> z
%

% Using function of RobustKnock
% [A_w,Ay_w ,B_w,C_w, lb_w, ub_w, wSize, wZs]=seperateTransposeJoin1(A_pk4(:,1:end-1), A_pk4(:,end), b_pk4,c_p1 , 1e6, 1, zSize);
% [A_w,Ay_w ,B_w,C_w, lb_w, ub_w, wSize, wZs]=seperateTransposeJoin1(A_pk4(:,1:end-1), A_pk4(:,end), b_pk4,c_p1 , 1,     1,    30,       1e6,     1,           3);
% ctype_dk1 = [repmat('C',1,size(A_w,2)) repmat('B',1,size(Ay_w,2)) ];
% [fv1,opt_ub_dk1,~] = cplexmilp([C_w;0],[A_w,Ay_w],B_w,[],[],[],[],[],lb_w,ub_w,ctype_dk1);
% 
[A_w,Ay_w ,B_w,C_w, lb_w, ub_w, wSize, wZs]=seperateTransposeJoin1(A_pk3(:,1:end-1), A_pk3(:,end), b_pk3,c_p1 ,  1e6, 1, 1);
A_rk = [A_w,Ay_w];
b_rk = B_w;
c_rk = C_w;
lb_rk = lb_w;
ub_rk = ub_w;
ctype_rk = [repmat('C',1,size(A_w,2)),repmat('B',1,size(Ay_w,2))];

[~,opt_rk,~] = intlinprog([c_rk; 0], find(ctype_rk=='B'), A_rk,b_rk,[],[],lb_rk,ub_rk,[],intlinprog_options);

  
%% 6. Dual problem with knockable reaction 6 / MILP
% Should give a higher maximum production at maximum biomass when reaction is knocked out
A_dk1 = -A_pk2';
b_dk1 = c_pk2;
c_dk1 = b_pk2;
lb_dk1 = zeros(size(A_dk1,2),1);
ub_dk1 = inf(size(A_dk1,2),1);
ctype_dk1 = repmat('C',1,size(A_dk1,2));

% Ideen hier
