%%  Basic RBC model %% 
%%  By Ryo Kato, modified January, 2004 
%%  variable elastic labor supply

%%  This code is written by Ryo Kato (kato.13@osu.edu)
%%  You can use/rewrite this code without my permission.
%%  I will not be responsible for any damage/liability 
%%  induced by running this conde.

clear all;

%%  ATTENTION!!
%%  Labor disutility ==>   -myu*(lt^(1+lamda))/(1+lamda)


%% ------------------- [1] Parameter proc  ------------------------
  sigma = 1.5;   % CRRA
  alpha = 0.3;   % Cob-Dag
  myu = 1;
  beta = 0.99;   
  delta = 0.025; 
  lamda = 2;   % labor supply elasticity >1
  phi = 0.8;     % AR(1) in tech


%% --------------------- [2] Steady State proc >> -----------------------
 % SS capital & ss labor
 % (1) real rate  (By SS euler)
        kls = (((1/beta)+delta-1)/alpha)^(1/(alpha-1));
 % (2) wage
        wstar = (1-alpha)*(kls)^alpha;      
 % (3) Labor and goods market clear
        clstar = kls^alpha - delta*kls;
        lstar = ((wstar/myu)*(clstar^(-sigma)))^(1/(lamda+sigma)) ;
        kstar = kls*lstar ;     
        cstar = clstar*lstar ;
        vstar = 1;
        Ystar = (kstar^alpha)*(lstar^(1-alpha));
     
    ssCKoLY = [ cstar  kstar ; lstar  Ystar ]  % show SS values
   
  
%% --------------------------[2] MODEL proc-----------------------------%%  
   % Define endogenous vars  ('a' denotes t+1 values)
       syms  ct   kt  lt  rt  vt  ca   ka    la  ra  va  ; 
   % Eliminate Price  
       ra = (va*alpha*(ka/la)^(alpha-1)) ;
       wt = (1-alpha)*vt*(kt/lt)^alpha;
   % Optimal Conditions  & state transition
       labor   = lt^lamda-wt/(myu*ct^sigma) ;       % LS = LD 
       euler   = ct^(-sigma) -(ca^(-sigma))*beta*(1+ra-delta) 	  ;     % C-Euler
       capital = ka - (1-delta)*kt-vt*(kt^alpha)*(lt^(1-alpha))+ct ; % K-trans
       tech    = va - phi*vt;        
       
       optcon  = [labor ; euler  ; capital  ; tech];    
    
    % GDP (Optional) %
         Yt  = vt*(kt^alpha)*(lt^(1-alpha));

       
%% ------------------  [3] Linearization proc  ------------------------%%     
  % Differentiation 
     xx = [la  ca  ka  va  lt  ct  kt   vt] ; % define vars
     jopt = jacobian(optcon,xx);
     
     % For GDP (Optipnal)    
        xr = [lt  kt  vt]
        jy = jacobian(Yt,xr);

   % Evaluate each derivative
     kt = kstar;      ka=kstar;
     ct = cstar;     ca = cstar;
     la = lstar;     lt = lstar;
     vt = vstar;   va = vstar;
     
   % Define Linear Coefficients  
     coef = eval(jopt);
     
     % For GDP (optional) 
        coy  = eval(jy);
   % In terms of % deviations from ss    
    vo = [ lstar   cstar   kstar  vstar  ] ;
    TW = [ vo ; vo  ; vo ; vo ] ;
 
    B =  [ -coef(:,1:4)  ].*TW ;
    C =  [  coef(:,5:8)  ].*TW ;
 % B[c(t+1)  l(t+1)  k(t+1)  z(t+1)] = C[c(t)  l(t)  k(t)  z(t)]
    A = inv(C)*B  %(Linearized reduced form )
   
   % For GDP( optional) 
      ve  = [lstar  kstar  vstar  ];
      NOM = [Ystar  Ystar  Ystar  ];
      PPX = coy.*ve./NOM;

%% =========== [4] Solution proc (Do NOT touch!!) ============== %%
%  EIGEN DECOMPOSITION
   [W V] = eig(A);
   Q = inv(W);
   W*V*Q
   theta = diag(V)
  % Extract stable vectors
      SQ = [];    jw = 1;
      for j = 1:length(theta)
      	if abs(theta(j)) >1.000000001
         	SQ(jw,:) = Q(j,:);
         	jw = jw+1;
      	end
      end
   % Extract unstable vectors
     UQ = [];   jjw = 1;
      for jj = 1:length(theta)
      	if abs(theta(jj))<0.9999999999
         	UQ(jjw,:) = Q(jj,:);
         	jjw = jjw+1;
      	end
      end
    % Extract stable roots 
      VLL = [];    jjjw = 1;
      for jjj = 1:length(theta)
      	if abs(theta(jjj)) >1.0000000001
            VLL(jjjw,:) = theta(jjj,:); 
            jjjw = jjjw+1;
      	end
      end
      
   % Show Eigen Vectors on U-S Roots  
     UQ;  % n x n+k
     SQ;  % k x n+k 
     
  % [3] ELIMINATING UNSTABLE VECTORS
      k = min(size(SQ));   % # of predetermined vars
      n = min(size(UQ));   % # of jump vars
      nk = [n k]
   % Stable V (eig mat)
   VL = inv(diag(VLL))
   
   % Elements in Q
    PA = UQ(1:n,1:n);    PB = UQ(1:n,n+1:n+k);
    PC = SQ(1:k,1:n);    PD = SQ(1:k,n+1:n+k);
    P = -inv(PA)*PB ; % X(t) = P*S(t)
    PE = PC*P+PD ;
   
   % SOLUTION
    PX = inv(PE)*VL*PE;
    AA = real(PX)   

    
%% ------------------ [5]  SIMULATION proc  ----------------- %%    
% [4] TIME&INITIAL VALUES 

     t = 24;			% Time span   
    
   % Initial Values 
   % state var + e
   S1 = [ 0   0.06 ]' ;

    
% [5] SIMULATION 
      Ss = S1;
      S = zeros(t,k) ;  
        for i = 1:t
          q = AA*Ss  ;
   	      S(i,:) = q';
          Ss = S(i,:)';
     	end
      SY = [S1' ;S]  ; 
      X = (real(P)*SY')';
      
   % Re-definition   
     ci = X(:,1);
     li = X(:,2);
     ki = SY(:,1);
     vi = SY(:,2);
      
     XI = [li ki vi];
     Yi = (PPX*XI')';      
% [6] DRAWING FIGURES
   figure(1)
     subplot(2,1,1);
       plot( Yi );
       ylabel('Yt')
       xlabel('time');
     subplot(2, 1,2);
       plot([ci li] );
       ylabel('C(blue)  L(grn)')
       xlabel('time');
     







