close all; clear; clc;

% constants
unitContourSize = 1000;
gridSize = 0.05;
gridSize = 0.5;

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
fig1 = figure; hold on;
plot (A(:,1), A(:,2), 'r.');
plot (B(:,1), B(:,2), 'g.');
plot (A_unitContour(:,1), A_unitContour(:,2), 'r-');
plot (B_unitContour(:,1), B_unitContour(:,2), 'g-');
axis equal;

% Case 2 plot
fig2 = figure; hold on;
plot (C(:,1), C(:,2), 'r.');
plot (D(:,1), D(:,2), 'g.');
plot (E(:,1), E(:,2), 'b.');
plot (C_unitContour(:,1), C_unitContour(:,2), 'r-');
plot (D_unitContour(:,1), D_unitContour(:,2), 'g-');
plot (E_unitContour(:,1), E_unitContour(:,2), 'b-');
axis equal;

%===============================================================================
% 3. CLASSIFIERS
%===============================================================================
% Grid prep
xMin1 = floor(min([A(:,1);B(:,1)])/5)*5;
xMax1 =  ceil(max([A(:,1);B(:,1)])/5)*5;
yMin1 = floor(min([A(:,2);B(:,2)])/5)*5;
yMax1 =  ceil(max([A(:,2);B(:,2)])/5)*5;
xAxis1 = xMin1:gridSize:xMax1;
yAxis1 = yMin1:gridSize:yMax1;
blankGrid1 = zeros(length(yAxis1),length(xAxis1));

xMin2 = floor(min([C(:,1);D(:,1);E(:,1)])/5)*5;
xMax2 =  ceil(max([C(:,1);D(:,1);E(:,1)])/5)*5;
yMin2 = floor(min([C(:,2);D(:,2);E(:,2)])/5)*5;
yMax2 =  ceil(max([C(:,2);D(:,2);E(:,2)])/5)*5;
xAxis2 = xMin2:gridSize:xMax2;
yAxis2 = yMin2:gridSize:yMax2;
blankGrid2 = zeros(length(yAxis2),length(xAxis2));

mu_A_index = [find(mu_A(2) == yAxis1) find(mu_A(1) == xAxis1)];
mu_B_index = [find(mu_B(2) == yAxis1) find(mu_B(1) == xAxis1)];
mu_C_index = [find(mu_C(2) == yAxis2) find(mu_C(1) == xAxis2)];
mu_D_index = [find(mu_D(2) == yAxis2) find(mu_D(1) == xAxis2)];
mu_E_index = [find(mu_E(2) == yAxis2) find(mu_E(1) == xAxis2)];

% MED
MEDGrid1 = blankGrid1;
for i = 1:size(MEDGrid1,1)
  for j = 1:size(MEDGrid1,2)
    [value, index] = min([-mu_A_index*[i j]'+0.5*mu_A_index*mu_A_index';
                          -mu_B_index*[i j]'+0.5*mu_B_index*mu_B_index']);
    MEDGrid1(i,j) = index-1;
  end
end

MEDGrid2 = blankGrid2;
for i = 1:size(MEDGrid2,1)
  for j = 1:size(MEDGrid2,2)
    [value, index] = min([-mu_C_index*[i j]'+0.5*mu_C_index*mu_C_index';
                          -mu_D_index*[i j]'+0.5*mu_D_index*mu_D_index';
                          -mu_E_index*[i j]'+0.5*mu_E_index*mu_E_index']);
    MEDGrid2(i,j) = index-1;
  end
end

figure(fig1);
contour(xAxis1,yAxis1,MEDGrid1,1, '-k');

figure(fig2);
contour(xAxis2,yAxis2,MEDGrid2,2, '-k');

% GED
figure(fig1);

figure(fig2);
% MAP
figure(fig1);

figure(fig2);
% NN
figure(fig1);

figure(fig2);

