clear all
%% 7. Bilevel mit Knockouts
%Anzahl für minimale und maximale
%Knockouts werden entgegengenommen
%ohne split
    %7.1 Primal with knockout
    load('duality_model.mat');
    cnap = model;
    idx_r6 = find(ismember(cellstr(cnap.reacID),'r_6'));
   
    cnap.reacMax(idx_r6) = 100; %r_6 ist i.A. eingeschalten

    A  = cnap.stoichMat;
    b  = zeros(cnap.nums,1);
    lb = cnap.reacMin;
    %lb(4)=0.5; %mindestwachstum
    ub = cnap.reacMax;
    c  = cnap.objFunc;
    n  = size(cnap.stoichMat,2);
    k_max=6; %Anzahl Knockouts
    k_min=0;
  
    
    %boundaries in A packen
    A_ub=eye(n); 
    A_lb=-eye(n);
    
    %Reshape into Ax<=b
    A_k2=[A; -A; A_ub; A_lb];

   
    lb1=max(zeros(n,1),lb); %vgl split
    lb2=max(zeros(n,1),-ub);
    ub1=max(zeros(n,1),ub);
    ub2=max(zeros(n,1),-lb);
    b_k=[b;-b; ub1; -lb1]; 
   
    %7.2 Dual with knockout

    c_d1=[b;-b; ub1; ub2; -lb1; -lb2]; %aus Splitting
    b_d1=[c; -c];
    A_plk=[A,-A]; 
    A_ub=eye(2*n); 
    A_lb=-eye(2*n);
    A_d1=-[A_plk; -A_plk; A_ub; A_lb]'; 

   
    M=10000;
    I=-M*eye(n);
    I2=[I;I];
    A_dK=[A_d1 I2]; %n spalten für n z Variablen hinzufügen
   

    c_dK=[c_d1; zeros(n,1)] ; %Zielfunktion um 0 für z erweitern
    b_dK=b_d1; %aus Dualem übernehmen

    
    %7.3 Bilevel erstellen
    
    Z=zeros(size(A_k2,1), size(A_d1,2));
    %(1-z)*ub => +z*Ub erstellen
    I3=diag(ub1);
   
    I5=-diag(lb1);
   
  
    I7=[zeros(2*size(A,1),n); I3; I5]; %(ub, lb, -lb, -ub), max{0,ub..} 

    
    
    A_ineqbi1_K=[A_k2, Z, I7 ; zeros(size(A_dK,1), size(A_k2,2)), A_dK; zeros(1, size(A_dK,2)) ones(1,n);zeros(1,size(A_dK,2)) -ones(1,n)]; %primal,dual Problem; zwei Zeilen für k_max, k_min)
        
    
    b_ineqbi1_K=[b_k;b_d1;k_max; -k_min;]; 
    
    
    z_fest_lb=zeros(n,1);
    %z_fest_lb(10)=1; %r_6 ausschalten
    z_fest_ub=ones(n,1); 
    %z_fest_ub(1)=0;%r_S darf nicht ausgeschalten werden
    %z_fest_ub(4)=0;%r_Bm darf nicht ausgeschalten werden
    %z_fest_ub(10)=0;
    
    %hier z festsetzen, z=0 acitve (7 BM), z=1 not active-mit knockout (3,5 BM)
    lb_K=[zeros(length(c_k),1) ;zeros(length(c_d1),1); z_fest_lb]; %-inf, damit primale Variablen frei (wegen splitting eigentlich nicht notwendig
    ub_K=[inf(length(c_k),1) ; inf(length(c_d1),1); z_fest_ub];

    %constraint to use strong duality
    Aeq_K=[c' c_dK']; %c aus 0. c_dK aus aus 6.strong duality
    beq_K=[0]; 

    c_prod_K=zeros(1, size (A_ineqbi1_K,2));
    c_prod_K(2)=-1; 
    c_prod_K(2+n)=1; 

    intcon=(size(A_ineqbi1_K,2)-(n-1):size(A_ineqbi1_K,2));
    
[R, valprod_K,~,~]=intlinprog(c_prod_K,intcon,A_ineqbi1_K, b_ineqbi1_K, Aeq_K, beq_K, lb_K, ub_K);

R;
