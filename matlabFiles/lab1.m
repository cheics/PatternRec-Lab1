close all; clear; clc;

% constants
unitContourSize = 1000;
gridSize = 0.05;
% gridSize = 0.5;

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
%% 2. GENERATING CLUSTERS
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

% ===============================================================================
%% 3. CLASSIFIERS
% ===============================================================================
% Grid prep
[xVals_AB, yVals_AB, MED_AB] = gridPrep(gridSize, A, B);
[xVals_CDE, yVals_CDE, MED_CDE] = gridPrep(gridSize, C, D, E);

xyGrid_AB=zeros(size(MED_AB,1), size(MED_AB,2), 2);
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


%% 3.1 MED Class
% % MED_AB
% for i = 1:size(MED_AB,1)
%   for j = 1:size(MED_AB,2)
% 		MED_AB(i,j)= MED_Class([xVals_AB(j), yVals_AB(i)], ...
% 			mu_A, mu_B);
%   end
% end
% 
% % MED_CDE
% for i = 1:size(MED_CDE,1)
%   for j = 1:size(MED_CDE,2)
% 		MED_CDE(i,j)= MED_Class([xVals_CDE(j), yVals_CDE(i)], ...
% 			mu_C, mu_D, mu_E);
%   end
% end

% figure(fig1);
% contour(xVals_AB,yVals_AB,MED_AB,1, '-k');
% figure(fig2);
% contour(xVals_CDE,yVals_CDE,MED_CDE,2, '-k');

%% 3.2 NN Class
% % NN_AB
% NN_class_AB=knnclassify(xyGrid_list_AB, ...
% 	vertcat(A, B), ...
% 	vertcat(ones(length(A),1), ones(length(B),1).*2) ...
% );
% NN_class_AB=reshape(NN_class_AB,size(MED_AB,1),size(MED_AB,2));
% 
% % NN_CDE
% NN_class_CDE=knnclassify(xyGrid_list_CDE, ...
% 	vertcat(C, D, E), ...
% 	vertcat(ones(length(C),1), ones(length(D),1).*2, ones(length(E),1).*3), ...
% );
% NN_class_CDE=reshape(NN_class_CDE,size(MED_CDE,1),size(MED_CDE,2));
% 
% figure(fig1);
% contour(xVals_AB,yVals_AB,NN_class_AB,1, '-k');
% figure(fig2);
% contour(xVals_CDE,yVals_CDE,NN_class_CDE,2, '-k');

%% 3.3 5NN Class
% % NN5_AB
% NN5_class_AB=knnclassify(xyGrid_list_AB, ...
% 	vertcat(A, B), ...
% 	vertcat(ones(length(A),1), ones(length(B),1).*2), ...
% 	5 ...
% );
% NN5_class_AB=reshape(NN5_class_AB,size(MED_AB,1),size(MED_AB,2));
% 
% % NN5_CDE
% NN5_class_CDE=knnclassify(xyGrid_list_CDE, ...
% 	vertcat(C, D, E), ...
% 	vertcat(ones(length(C),1), ones(length(D),1).*2, ones(length(E),1).*3), ...
% 	5 ...
% );
% NN5_class_CDE=reshape(NN5_class_CDE,size(MED_CDE,1),size(MED_CDE,2));
% 
% figure(fig1);
% contour(xVals_AB,yVals_AB,NN5_class_AB);
% figure(fig2);
% contour(xVals_CDE,yVals_CDE,NN5_class_CDE);


% 
% % GED
% figure(fig1);
% 
% figure(fig2);
% % MAP
% figure(fig1);
% 
% figure(fig2);

