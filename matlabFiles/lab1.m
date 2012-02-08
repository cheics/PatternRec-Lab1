close all; clear; clc;

VISIBLE_FLAG = 1;

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
