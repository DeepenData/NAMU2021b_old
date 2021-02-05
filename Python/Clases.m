



A[e] -> A[c]
A[c] + atp[c] -> B[c] + adp[c] + pi[c]
B[c] + 2 adp[c]+ pi[c] + 3 nad[c] + 3 h[c]  -> 2 atp[c] + 3 nadh[c] + C[c]
C[c] + 2 nad[c] + 2 h[c] -> 0.2 C[c] + 2 nadh[c] 
C[c] + adp[c] + pi[c] -> atp[c] + 3 D[c]
C[c] + 2 nadh[c] -> 3 E[c] + 2 nad[c] + 2 h[c]
nadh[c] + o2[c]+ 2 adp[c] + 2 pi[c] -> 2 atp[c] + nad[c] + h[c]
atp[c] -> adp[c] + pi[c] 
C[c] + 10 atp[c] -> biomass[c] + 10 adp[c] + 10 pi[c] 
o2[e] -> o2[c] 
C[e] <-> C[c]
D[e] <-> D[c]
E[e] <-> E[c]






%------------------------------
%cd('C:\Users\alejandro\Documents\MATLAB\cobra');
initCobraToolbox;
reactionFormulas = {...
'<=>  A[c]',...%A_uptake
'A[c] + atp[c] -> B[c]',...%r1
'B[c]  -> 2 atp[c] + 3 nadh[c] + C[c]',...%r2
'C[c] -> 2 nadh[c] + 0.2 C[c]',... %r3
'C[c]  -> atp[c] + 3 D[c]',... %r4
'C[c] + 2 nadh[c] -> 3 E[c]',... %r5
'nadh[c] + o2[c] -> 2 atp[c]',...%r_res
'atp[c] ->',... %r_m 
'C[c] + 10 atp[c] ->',... % rbio objetivo
' <=> o2[c]',... %o2_uptake
' <=> C[c]',...  %C_exchange
' <=> D[c]',...  %D_exchange
' <=> E[c]'};    %E_exchange
reactionNames = {'A_uptake', 'r1', 'r2', 'r3', 'r4','r5','r_res','r_m','rbio',...
    'o2_uptake','C_exchange','D_exchange','E_exchange'};
%'A_uptake','r1','r2','r3','r4','r5','r_res','r_m','rbio' 'o2_uptake','C_exchange','D_exchange','E_exchange'
lowerBounds = [-10, 0,  0,  0,  0,  0,  0,  5,  0, -100, -100, -100, -100];
upperBounds = [10, 20, 20, 20, 20, 20, 20, 20, 20,  100,  100,  100,  100];
Toy_Model = createModel(reactionNames,reactionNames,reactionFormulas,[],...
    lowerBounds,upperBounds)
Toy_Model = changeObjective(Toy_Model,'rbio');
save('C:\Users\alejandro\Dropbox\2017\Modelación y simulación de procesos\Models.mat','Toy_Model');
clear,clc
load('C:\Users\alejandro\Dropbox\2017\Modelación y simulación de procesos\Models.mat');


%c=Toy_Model.c;osense=-1;A=Toy_Model.S;b=Toy_Model.b;lb=Toy_Model.lb;ub=Toy_Model.ub;
% [x,f,origStat,output,lambda] = linprog(c*osense,[],[],A,b,lb,ub);

%INPUT
% LPproblem Structure containing the following fields describing the LP
% problem to be solved
%  A      LHS matrix, left-hand side matrix. Es la matriz estequiométrica
%  b      RHS vector, right-hand side vector.
%  c      Objective coeff vector
%  lb     Lower bound vector
%  ub     Upper bound vector
%  osense Objective sense (-1 max, +1 min)
%  csense Constraint senses, a string containting the constraint sense for
%         each row in A ('E', equality, 'G' greater than, 'L' less than).

Toy_Model

%Definir los límites para los consumos de sustrato con COBRA

Toy_Model = changeRxnBounds(Toy_Model,'A_uptake',0,'l'); %lower bound
Toy_Model = changeRxnBounds(Toy_Model,'A_uptake',10,'u');%upper bound
Toy_Model = changeRxnBounds(Toy_Model,'o2_uptake',0,'l');%lower bound
Toy_Model = changeRxnBounds(Toy_Model,'o2_uptake',10,'u');%upper bound
Toy_Model = changeRxnBounds(Toy_Model,'r_m',10,'u');%upper bound

Toy_Model = changeRxnBounds(Toy_Model,'C_exchange',100,'u');%upper bound
Toy_Model = changeRxnBounds(Toy_Model,'C_exchange',-100,'l');%lower bound

Toy_Model = changeRxnBounds(Toy_Model,'D_exchange',100,'u');%upper bound
Toy_Model = changeRxnBounds(Toy_Model,'D_exchange',0,'l');

Toy_Model = changeRxnBounds(Toy_Model,'E_exchange',100,'u');%upper bound
Toy_Model = changeRxnBounds(Toy_Model,'E_exchange',-100,'l');

opt = optimizeCbModel(Toy_Model,'max') % f: 4.9000


r=0;clear growth_rate fluxes  %contador en cero
for i=0:0.5:20; %rango de valores para el consumo de carbono
    r=r+1;   %contador 
    z=0;     %segundo contador en cero
    for j=0:0.5:10; % rango de valores para el consumo de oxígeno
        z=z+1;   %segundo contador
  Toy_Model = changeRxnBounds(Toy_Model,'A_uptake',i,'u');%upper bound
  Toy_Model = changeRxnBounds(Toy_Model,'A_uptake',0,'l');%lower bound
  Toy_Model = changeRxnBounds(Toy_Model,'o2_uptake',j,'u');%upper bound
  Toy_Model = changeRxnBounds(Toy_Model,'o2_uptake',0,'l');%lower bound
  try
  opt = optimizeCbModel(Toy_Model,'max'); %Optimización para resolver
  growth_rate(r,z)= opt.f ; %Valor de la función objetivo
  fluxes(r,z,:)=opt.x; %Valor de todos los flujos
  size(fluxes)
  end
    end
end
%             J(o2_uptake),I(A_uptake),var
figure(1);box on;
h(1)=subplot(2,2,1);surf(0:0.5:10,0:0.5:20,growth_rate,'EdgeColor','none','FaceLighting','phong'); xlabel('o2_uptake');ylabel('A_uptake');zlabel('growth_rate');

h(2)=subplot(2,2,2);surf(0:0.5:10,0:0.5:20,fluxes(:,:,find(ismember(Toy_Model.rxns,'C_exchange'))...
),'EdgeColor','none','FaceLighting','phong'); xlabel('o2_uptake');ylabel('A_uptake');zlabel('product');

h(3)=subplot(2,2,3);surf(0:0.5:10,0:0.5:20,fluxes(:,:,find(ismember(Toy_Model.rxns,'D_exchange'))...
),'EdgeColor','none','FaceLighting','phong'); xlabel('o2_uptake');ylabel('A_uptake');


h(4)=subplot(2,2,4);surf(0:0.5:10,0:0.5:20,fluxes(:,:,find(ismember(Toy_Model.rxns,'E_exchange'))...
),'EdgeColor','none','FaceLighting','phong'); xlabel('o2_uptake');ylabel('A_uptake');



%

reactionFormulas = {...
'->  A[c]',...%A_uptake
'A[c] + atp[c] -> B[c]',...%r1
'B[c]  -> 2 atp[c] + 3 nadh[c] + C[c]',...%r2
'C[c] -> 2 nadh[c] + 0.2 C[c]',... %r3
'C[c]  -> atp[c] + 3 D[c]',... %r4
'C[c] + 2 nadh[c] -> 3 E[c]',... %r5
'nadh[c] + o2[c] -> 2 atp[c]',...%r_res
'atp[c] ->',... %r_m
'C[c] + 10 atp[c] ->',... % rbio objetivo
' -> o2[c]',... %o2_uptake
'C[c] <=>',...  %C_exchange
'D[c] <=>',...  %D_exchange
'E[c] <=>'};    %E_exchange
reactionNames = {'A_uptake', 'r1', 'r2', 'r3', 'r4','r5','r_res','r_m','rbio',...
    'o2_uptake','C_exchange','D_exchange','E_exchange'};
lowerBounds = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -20, -20, -20];
upperBounds = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20];
Toy_Model2 = createModel(reactionNames,reactionNames,reactionFormulas,[],...
    lowerBounds,upperBounds);Toy_Model2.mets;Toy_Model2.rxns;
Toy_Model2 = changeObjective(Toy_Model2,'rbio');
Toy_Model2 = changeRxnBounds(Toy_Model2,'A_uptake',1,'l'); 
Toy_Model2 = changeRxnBounds(Toy_Model2,'A_uptake',10,'u'); 
Toy_Model2 = changeRxnBounds(Toy_Model2,'o2_uptake',0,'l'); 
Toy_Model2 = changeRxnBounds(Toy_Model2,'o2_uptake',10,'u'); %-5.85
opt2 = optimizeCbModel(Toy_Model2,'max')

[opt2.x([1]);Toy_Model2.rxns([1])]




%Flux balance analysis (FBA)
% Linear Programming (LP) problem of the form: 
%  c: vector columna             
%                     min c'*v     Para el caso del toy_model c' sólo selecciona a rbio      
%c: coeficientes estequiométricos de la función objetivo (biomasa). Tiene la longitud
%del número de flujos, selecciona los flujos que forman parte de la función
%objetivo.
%': transposición 
%v: el vector de flujos o reacciones. 
%Para maximizar debo poner un signo negativo, es decir, cambiar el sentido: min -1*c'*v 

%              
%    subject to S*v = b         restricciones de igualdad
%              lb <= v <= ub     restricciones de desigualdad, define la
%              reversibilidad, por ejemplo: irreversible lb_i=0, ub=10

%
%
%INPUT
% LPproblem Structure containing the following fields describing the LP
% problem to be solved
%  A      LHS matrix
%  b      RHS vector
%  c      Objective coeff vector
%  lb     Lower bound vector
%  ub     Upper bound vector
%  osense Objective sense (-1 max, +1 min)
%  csense Constraint senses, a string containting the constraint sense for
%         each row in A ('E', equality, 'G' greater than, 'L' less than).






c=Toy_Model2.c;
osense=-1;
A=Toy_Model2.S;b=Toy_Model2.b;lb=Toy_Model2.lb;ub=Toy_Model2.ub

 [x,f,origStat,output,lambda] = linprog(c*osense,[],[],A,b,lb,ub)






['Phase',model.rxns([573 518 1010 197 440])';{'Control',sol_1.x(573),sol_1.x(518)...
    ,sol_1.x(1010),sol_1.f,sol_1.x(440)};{'Inhibi',sol_2.x(573),sol_2.x(518)...
    ,sol_2.x(1010),sol_2.f,sol_2.x(440)};{'Difference',dif_inh(573),dif_inh(518)...
    ,dif_inh(1010),dif_inh(197),dif_inh(440)};...
    {'alt_shift',sol_s.x(573),sol_s.x(518)...
    ,sol_s.x(1010),sol_s.x(197),sol_s.x(440)}]'


%----------------------------------------------------------------------------

% To create this model, we have to define the reactions:
reactionFormulas = {'glc-D[e]  -> glc-D',...
    'glc-D + atp  -> H + adp + g6p',...
    'g6p  <=> f6p',...
    'atp + f6p  -> H + adp + fdp',...
    'fdp + h2o  -> f6p + pi',...
    'fdp  -> g3p + dhap',...
    'dhap  -> g3p'};
reactionNames = {'GLCt1', 'HEX1', 'PGI', 'PFK', 'FBP', 'FBA', 'TPI'};
lowerBounds = [-20, 0, -20, 0, 0, -20, -20];
upperBounds = [20, 20, 20, 20, 20, 20, 20];

%2017
%glycolysisModel = createModel(reactionNames, reactionNames, reactionFormulas,...
%                             'lowerBoundList', lowerBounds, 'upperBoundList', upperBounds);

%2012
glycolysisModel = createModel(reactionNames,reactionNames,reactionFormulas,[],...
    lowerBounds,upperBounds);


load('C:\Users\alejandro\Dropbox\2017\Modelación y simulación de procesos\iRC1080.mat');%savepath();
global iRC1080
printRxnFormula(iRC1080,iRC1080.rxns(find(iRC1080.c)));%shown objective function 
printRxnFormula(iRC1080,iRC1080.rxns(find(~cellfun(@isempty,regexpi(iRC1080.rxnNames,'biomass')))));%shown biomass rxns
iRC1080.c(find(ismember(iRC1080.rxns,'BIOMASS_Chlamy_hetero')))=1;%define objective 
iRC1080.c(find(ismember(iRC1080.rxns,'BIOMASS_Chlamy_mixo')))=1;%define objective 
iRC1080.c(find(ismember(iRC1080.rxns,'BIOMASS_Chlamy_auto')))=1;%define objective 

iRC1080=changeRxnBounds(iRC1080,{'EX_co2_e';'EX_o2_e';'EX_photonVis_e'},[-1.1 -1.1 -2000],'l');
iRC1080=changeRxnBounds(iRC1080,{'EX_co2_e';'EX_o2_e';'EX_photonVis_e'},[1 100 0],'u');
opt = optimizeCbModel(iRC1080,'max')
table(...
[printRxnFormula(iRC1080,iRC1080.rxns(exc))],[iRC1080.lb(exc)],[iRC1080.ub(exc)],opt.x(exc))




[Exc1, Upt1] = findExcRxns(iRC1080); %find exchange and uptake reactions
exc=union(find(Exc1),find(Upt1));%get all exchange and uptake reactions together
printRxnFormula(iRC1080,iRC1080.rxns(exc));%shown all exc
exc_rxns=iRC1080.rxns(exc);
printRxnFormula(iRC1080,iRC1080.rxns(find(~cellfun(@isempty,regexpi(iRC1080.rxns,'starch')))));%search among all rxns
printRxnFormula(iRC1080,exc_rxns(find(~cellfun(@isempty,regexpi(exc_rxns,'o2_e')))));%search inside exc_rxns

table(...
[printRxnFormula(iRC1080,iRC1080.rxns(exc))],[iRC1080.lb(exc)],[iRC1080.ub(exc)],opt.x(exc))
r=0;clear growth_rate fluxes  %contador en cero
for i=0:0.005:0.1; %rango de valores para el consumo de carbono
    r=r+1;   %contador 
    z=0;     %segundo contador en cero
    for j=0:0.005:0.1; % rango de valores para el consumo de oxígeno
        z=z+1;   %segundo contador
  iRC1080 = changeRxnBounds(iRC1080,'EX_co2_e',i,'u');%upper bound
  iRC1080 = changeRxnBounds(iRC1080,'EX_co2_e',i,'l');%lower bound
  iRC1080 = changeRxnBounds(iRC1080,'EX_o2_e',-j,'u');%upper bound
  iRC1080 = changeRxnBounds(iRC1080,'EX_o2_e',-j,'l');%lower bound
  try
  opt = optimizeCbModel(iRC1080,'max'); %Optimización para resolver
  growth_rate(r,z)= opt.f ; %Valor de la función objetivo
  fluxes(r,z,:)=opt.x; %Valor de todos los flujos
  size(fluxes)
  end
    end
end
%             J(o2_uptake),I(A_uptake),var
figure(1);box on;
h(1)=subplot(2,2,1);surf(0:0.005:0.1,0:0.005:0.1,growth_rate,'EdgeColor','none','FaceLighting','phong'); xlabel('o2_uptake');ylabel('A_uptake');zlabel('growth_rate');

h(2)=subplot(2,2,2);surf(0:0.5:10,0:0.5:20,fluxes(:,:,find(ismember(Toy_Model.rxns,'C_exchange'))...
),'EdgeColor','none','FaceLighting','phong'); xlabel('o2_uptake');ylabel('A_uptake');zlabel('product');

h(3)=subplot(2,2,3);surf(0:0.5:10,0:0.5:20,fluxes(:,:,find(ismember(Toy_Model.rxns,'D_exchange'))...
),'EdgeColor','none','FaceLighting','phong'); xlabel('o2_uptake');ylabel('A_uptake');


h(4)=subplot(2,2,4);surf(0:0.5:10,0:0.5:20,fluxes(:,:,find(ismember(Toy_Model.rxns,'E_exchange'))...
),'EdgeColor','none','FaceLighting','phong'); xlabel('o2_uptake');ylabel('A_uptake');








