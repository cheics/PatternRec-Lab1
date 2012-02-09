close all; clear; clc;

VISIBLE_FLAG = 0;

% constants
unitContourSize = 10000;
%gridSize = 0.02;
gridSize = 0.05;
% gridSize = 0.5;
fontSize = 10;
figSize = [0 0 6 4];% [something, something, width, height]
shadeVal=0.07;

% Cluster data
N_A = 200;
N_B = 200;
N_C = 100;
N_D = 200;
N_E = 150;

mu_A = [ 5 10];
mu_B = [10 15];
mu_C = [ 5 10];
mu_D = [15 10];
mu_E = [10  5];

Sigma_A = [ 8  0;  0  4];
Sigma_B = [ 8  0;  0  4];
Sigma_C = [ 8  4;  4 40];
Sigma_D = [ 8  0;  0  8];
Sigma_E = [10 -5; -5 20];

%===============================================================================
% 2. GENERATING CLUSTERS
%===============================================================================

% create clusters
A = gaussTransform(randn(N_A,2),mu_A,Sigma_A);
B = gaussTransform(randn(N_B,2),mu_B,Sigma_B);
C = gaussTransform(randn(N_C,2),mu_C,Sigma_C);
D = gaussTransform(randn(N_D,2),mu_D,Sigma_D);
E = gaussTransform(randn(N_E,2),mu_E,Sigma_E);

% create generic unit contour (i.e. a circle)
unitContour = [cos((1:unitContourSize)/unitContourSize*2*pi());sin((1:unitContourSize)/unitContourSize*2*pi())]';

% transform unit contour for each class
A_unitContour = gaussTransform(unitContour,mu_A,Sigma_A);
B_unitContour = gaussTransform(unitContour,mu_B,Sigma_B);
C_unitContour = gaussTransform(unitContour,mu_C,Sigma_C);
D_unitContour = gaussTransform(unitContour,mu_D,Sigma_D);
E_unitContour = gaussTransform(unitContour,mu_E,Sigma_E);

% Case 1 plot
if VISIBLE_FLAG ~= 0
  fig1 = figure('paperposition', figSize); hold on;
else
  fig1 = figure('visible','off','paperposition', figSize); hold on;
end
plot (A(:,1), A(:,2), 'r.');
plot (B(:,1), B(:,2), 'g.');
plot (A_unitContour(:,1), A_unitContour(:,2), 'r-');
plot (B_unitContour(:,1), B_unitContour(:,2), 'g-');
axis equal;
xlabel('Feature 1', 'fontsize', fontSize);
ylabel('Feature 2', 'fontsize', fontSize);
print -dpng -r300 'fig1a-AB_cluster'

% Case 2 plot
if VISIBLE_FLAG ~= 0
  fig2 = figure('paperposition', figSize); hold on;
else
  fig2 = figure('visible','off','paperposition', figSize); hold on;
end
plot (C(:,1), C(:,2), 'r.');
plot (D(:,1), D(:,2), 'g.');
plot (E(:,1), E(:,2), 'b.');
plot (C_unitContour(:,1), C_unitContour(:,2), 'r-');
plot (D_unitContour(:,1), D_unitContour(:,2), 'g-');
plot (E_unitContour(:,1), E_unitContour(:,2), 'b-');
axis equal;
xlabel('Feature 1', 'fontsize', fontSize);
ylabel('Feature 2', 'fontsize', fontSize);
print -dpng -r300 'fig1b-CDE_cluster'

%===============================================================================
% 3. CLASSIFIERS
%===============================================================================
% Grid prep
[xVals_AB, yVals_AB, MED_AB] = gridPrep(gridSize, A, B);
[xVals_CDE, yVals_CDE, MED_CDE] = gridPrep(gridSize, C, D, E);

xyGrid_AB=zeros(size(MED_AB,1), size(MED_AB,2), 2);

GED_AB=MED_AB;
GED_CDE=MED_CDE;
MAP_AB=MED_AB;
MAP_CDE=MED_CDE;

for j = 1:size(MED_AB,1)
	for i = 1:size(MED_AB,2)
		xyGrid_AB(j, i, :) = [xVals_AB(i), yVals_AB(j)];
	end
end
xyGrid_list_AB=reshape(xyGrid_AB,size(MED_AB,1)*size(MED_AB,2),2);

xyGrid_CDE=zeros(size(MED_CDE,1), size(MED_CDE,2), 2);
for j = 1:size(MED_CDE,1)
	for i = 1:size(MED_CDE,2)
		xyGrid_CDE(j, i, :) = [xVals_CDE(i), yVals_CDE(j)];
	end
end
xyGrid_list_CDE=reshape(xyGrid_CDE,size(MED_CDE,1)*size(MED_CDE,2),2);

%-------------------------------------------------------------------------------
% 3.1 MED Class
%-------------------------------------------------------------------------------
% MED_AB
for i = 1:size(MED_AB,1)
  for j = 1:size(MED_AB,2)
		MED_AB(i,j)= MED_Class([xVals_AB(j), yVals_AB(i)], ...
			mu_A, mu_B);
  end
end

% MED_CDE
for i = 1:size(MED_CDE,1)
  for j = 1:size(MED_CDE,2)
		MED_CDE(i,j)= MED_Class([xVals_CDE(j), yVals_CDE(i)], ...
			mu_C, mu_D, mu_E);
  end
end

figure(fig1);
if VISIBLE_FLAG == 0
  set(fig1,'visible','off');
end
[c,h] =contour(xVals_AB,yVals_AB,MED_AB,1,'-c');
%ch = get(h,'child'); alpha(ch,shadeVal);

figure(fig2);
if VISIBLE_FLAG == 0
  set(fig2,'visible','off');
end
[c,h] =contour(xVals_CDE,yVals_CDE,MED_CDE,2,'-c');
%ch = get(h,'child'); alpha(ch,shadeVal);

%-------------------------------------------------------------------------------
% 3.2 GED Class
%-------------------------------------------------------------------------------
% GED_AB
for i = 1:size(GED_AB,1)
  for j = 1:size(GED_AB,2)
		GED_AB(i,j)=GED_Class2([xVals_AB(j), yVals_AB(i)],mu_A,Sigma_A,...
            N_A,mu_B,Sigma_B,N_B);
  end
end

% GED_CDE
for i = 1:size(GED_CDE,1)
  for j = 1:size(GED_CDE,2)
		GED_CDE(i,j)=GED_Class3([xVals_CDE(j), yVals_CDE(i)],...
            mu_C,Sigma_C,N_C,mu_D,Sigma_D,N_D,mu_E,Sigma_E,N_E);
  end
end
 
figure(fig1);
if VISIBLE_FLAG == 0
  set(fig1,'visible','off');
end
[c,h] =contour(xVals_AB,yVals_AB,GED_AB,1, '-m');
%ch = get(h,'child'); alpha(ch,shadeVal);
 
figure(fig2);
if VISIBLE_FLAG == 0
  set(fig2,'visible','off');
end
[c,h] =contour(xVals_CDE,yVals_CDE,GED_CDE,2, '-m');
%ch = get(h,'child'); alpha(ch,shadeVal);

%-------------------------------------------------------------------------------
% 3.3 MAP Class
%-------------------------------------------------------------------------------
% MAP_AB
for i = 1:size(MAP_AB,1)
  for j = 1:size(MAP_AB,2)
	MAP_AB(i,j)=MAP_class2([xVals_AB(j), yVals_AB(i)], ...
		mu_A,Sigma_A,N_A,mu_B,Sigma_B,N_B);
  end
end
figure(fig1);
if VISIBLE_FLAG == 0
  set(fig1,'visible','off');
end
[c,h] =contour(xVals_AB,yVals_AB,MAP_AB,1,'-b');
%ch = get(h,'child'); alpha(ch,shadeVal);
print -dpng -r300 'fig2a-AB_MED_MICD_MAP'

% MAP_CDE
for i = 1:size(MAP_CDE,1)
  for j = 1:size(MAP_CDE,2)
	MAP_CDE(i,j)=MAP_class3([xVals_CDE(j), yVals_CDE(i)], ...
		mu_C,Sigma_C,N_C,mu_D,Sigma_D,N_D,mu_E,Sigma_E,N_E);
  end
end
figure(fig2);
if VISIBLE_FLAG == 0
  set(fig2,'visible','off');
end
[c,h] =contour(xVals_CDE,yVals_CDE,MAP_CDE,2,'-b');
%ch = get(h,'child'); alpha(ch,shadeVal);
print -dpng -r300 'fig2b-AB_MED_MICD_MAP'

%-------------------------------------------------------------------------------
% 3.4 NN Class
%-------------------------------------------------------------------------------
% NN_AB
NN_class_AB=knnclassify(xyGrid_list_AB, ...
	vertcat(A, B), ...
	vertcat(ones(length(A),1), ones(length(B),1).*2) ...
);
NN_class_AB=reshape(NN_class_AB,size(MED_AB,1),size(MED_AB,2));

% NN_CDE
NN_class_CDE=knnclassify(xyGrid_list_CDE, ...
	vertcat(C, D, E), ...
	vertcat(ones(length(C),1), ones(length(D),1).*2, ones(length(E),1).*3) ...
	);
NN_class_CDE=reshape(NN_class_CDE,size(MED_CDE,1),size(MED_CDE,2));
 
if VISIBLE_FLAG ~= 0
  fig1 = figure('paperposition', figSize); hold on;
else
  fig1 = figure('visible','off','paperposition', figSize); hold on;
end
plot (A(:,1), A(:,2), 'r.');
plot (B(:,1), B(:,2), 'g.');
plot (A_unitContour(:,1), A_unitContour(:,2), 'r-');
plot (B_unitContour(:,1), B_unitContour(:,2), 'g-');
axis equal;
xlabel('Feature 1', 'fontsize', fontSize);
ylabel('Feature 2', 'fontsize', fontSize);
[c,h]= contourf(xVals_AB,yVals_AB,NN_class_AB,1);
ch = get(h,'child'); alpha(ch,shadeVal);
print -dpng -r300 'fig3a-AB_NN'

if VISIBLE_FLAG ~= 0
  fig2 = figure('paperposition', figSize); hold on;
else
  fig2 = figure('visible','off','paperposition', figSize); hold on;
end
plot (C(:,1), C(:,2), 'r.');
plot (D(:,1), D(:,2), 'g.');
plot (E(:,1), E(:,2), 'b.');
plot (C_unitContour(:,1), C_unitContour(:,2), 'r-');
plot (D_unitContour(:,1), D_unitContour(:,2), 'g-');
plot (E_unitContour(:,1), E_unitContour(:,2), 'b-');
axis equal;
xlabel('Feature 1', 'fontsize', fontSize);
ylabel('Feature 2', 'fontsize', fontSize);
[c,h]=contourf(xVals_CDE,yVals_CDE,NN_class_CDE,2);
ch = get(h,'child'); alpha(ch,shadeVal);
print -dpng -r300 'fig3b-CDE_NN'

%-------------------------------------------------------------------------------
% 3.5 5NN Class
%-------------------------------------------------------------------------------
% % NN5_AB
NN5_class_AB=knnclassify(xyGrid_list_AB, ...
	vertcat(A, B), ...
	vertcat(ones(length(A),1), ones(length(B),1).*2), ...
	5 ...
);
NN5_class_AB=reshape(NN5_class_AB,size(MED_AB,1),size(MED_AB,2));

% NN5_CDE
NN5_class_CDE=knnclassify(xyGrid_list_CDE, ...
	vertcat(C, D, E), ...
	vertcat(ones(length(C),1), ones(length(D),1).*2, ones(length(E),1).*3), ...
	5 ...
);
NN5_class_CDE=reshape(NN5_class_CDE,size(MED_CDE,1),size(MED_CDE,2));

if VISIBLE_FLAG ~= 0
  fig1 = figure('paperposition', figSize); hold on;
else
  fig1 = figure('visible','off','paperposition', figSize); hold on;
end
plot (A(:,1), A(:,2), 'r.');
plot (B(:,1), B(:,2), 'g.');
plot (A_unitContour(:,1), A_unitContour(:,2), 'r-');
plot (B_unitContour(:,1), B_unitContour(:,2), 'g-');
axis equal;
xlabel('Feature 1', 'fontsize', fontSize);
ylabel('Feature 2', 'fontsize', fontSize);
[c,h] =contourf(xVals_AB,yVals_AB,NN5_class_AB,1);
ch = get(h,'child'); alpha(ch,shadeVal);
print -dpng -r300 'fig4a-AB_5NN'

if VISIBLE_FLAG ~= 0
  fig2 = figure('paperposition', figSize); hold on;
else
  fig2 = figure('visible','off','paperposition', figSize); hold on;
end
plot (C(:,1), C(:,2), 'r.');
plot (D(:,1), D(:,2), 'g.');
plot (E(:,1), E(:,2), 'b.');
plot (C_unitContour(:,1), C_unitContour(:,2), 'r-');
plot (D_unitContour(:,1), D_unitContour(:,2), 'g-');
plot (E_unitContour(:,1), E_unitContour(:,2), 'b-');
axis equal;
xlabel('Feature 1', 'fontsize', fontSize);
ylabel('Feature 2', 'fontsize', fontSize);
[c,h] =contourf(xVals_CDE,yVals_CDE,NN5_class_CDE,2);
ch = get(h,'child'); alpha(ch,shadeVal);
print -dpng -r300 'fig4b-CDE_5NN'

%===============================================================================
% 4. ERROR ANALYSIS
%===============================================================================
% Finding P(E)
% The following method focuses on implementing the discrete implementation of
% equation 4.24 on page 68 of the SYDE 372 course notes.
%-------------------------------------------------------------------------------

% Define binary arrays of class regions R_A, etc.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
R_A = zeros(size(MAP_AB));
R_B = R_A;
R_C = zeros(size(MAP_CDE));
R_D = R_C;
R_E = R_C;

R_A(find(MAP_AB == 1)) = 1;
R_B(find(MAP_AB == 2)) = 1;
R_C(find(MAP_CDE == 1)) = 1;
R_D(find(MAP_CDE == 2)) = 1;
R_E(find(MAP_CDE == 3)) = 1;

% Define P(A), etc.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
P_of_A = N_A/(N_A+N_B);
P_of_B = N_B/(N_A+N_B);
P_of_C = N_C/(N_C+N_D+N_E);
P_of_D = N_D/(N_C+N_D+N_E);
P_of_E = N_E/(N_C+N_D+N_E);

% Define P(x|A), etc. for each discrete grid point
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
P_of_xgivenA = zeros(size(MAP_AB));
P_of_xgivenB = P_of_xgivenA;
P_of_xgivenC = zeros(size(MAP_CDE));
P_of_xgivenD = P_of_xgivenC;
P_of_xgivenE = P_of_xgivenC;

% see SYDE 372 course notes, page 21, equation 2.32 for definition of
% multivariate Gaussian distribution
for j = 1:size(P_of_xgivenA,1)
	for i = 1:size(P_of_xgivenA,2)
		P_of_xgivenA(j, i) = exp(-0.5*([xVals_AB(i), yVals_AB(j)]-mu_A)*inv(Sigma_A)*([xVals_AB(i), yVals_AB(j)]-mu_A)');
	end
end
P_of_xgivenA = P_of_xgivenA/(2*pi()*sqrt(det(Sigma_A)));

for j = 1:size(P_of_xgivenB,1)
	for i = 1:size(P_of_xgivenB,2)
		P_of_xgivenB(j, i) = exp(-0.5*([xVals_AB(i), yVals_AB(j)]-mu_B)*inv(Sigma_B)*([xVals_AB(i), yVals_AB(j)]-mu_B)');
	end
end
P_of_xgivenB = P_of_xgivenB/(2*pi()*sqrt(det(Sigma_B)));

for j = 1:size(P_of_xgivenC,1)
	for i = 1:size(P_of_xgivenC,2)
		P_of_xgivenC(j, i) = exp(-0.5*([xVals_CDE(i), yVals_CDE(j)]-mu_C)*inv(Sigma_C)*([xVals_CDE(i), yVals_CDE(j)]-mu_C)');
	end
end
P_of_xgivenC = P_of_xgivenC/(2*pi()*sqrt(det(Sigma_C)));

for j = 1:size(P_of_xgivenD,1)
	for i = 1:size(P_of_xgivenD,2)
		P_of_xgivenD(j, i) = exp(-0.5*([xVals_CDE(i), yVals_CDE(j)]-mu_D)*inv(Sigma_D)*([xVals_CDE(i), yVals_CDE(j)]-mu_D)');
	end
end
P_of_xgivenD = P_of_xgivenD/(2*pi()*sqrt(det(Sigma_D)));

for j = 1:size(P_of_xgivenE,1)
	for i = 1:size(P_of_xgivenE,2)
		P_of_xgivenE(j, i) = exp(-0.5*([xVals_CDE(i), yVals_CDE(j)]-mu_E)*inv(Sigma_E)*([xVals_CDE(i), yVals_CDE(j)]-mu_E)');
	end
end
P_of_xgivenE = P_of_xgivenE/(2*pi()*sqrt(det(Sigma_E)));

% Find P(E|x) (i.e. P(E) before integrating/summing; see equation 4.20) for each
% grid point
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
P_of_Egivenx_AB= P_of_xgivenA*P_of_A.*R_B + P_of_xgivenB*P_of_B.*R_A;
P_of_Egivenx_CDE= (P_of_xgivenD*P_of_D+P_of_xgivenE*P_of_E).*R_C + (P_of_xgivenC*P_of_C+P_of_xgivenE*P_of_E).*R_D + (P_of_xgivenC*P_of_C+P_of_xgivenD*P_of_D).*R_E;

% Find P(E) by summing P(E|x)*gridSize^2
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
P_of_E_AB = sum(P_of_Egivenx_AB(:))*gridSize^2;
P_of_E_CDE = sum(P_of_Egivenx_CDE(:))*gridSize^2;
