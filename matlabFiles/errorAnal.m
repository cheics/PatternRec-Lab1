dataSet_AB=vertcat([A, ones(length(A),1)*1], [B, ones(length(B),1)*2]);
dataSet_CDE=vertcat([C, ones(length(C),1)*1], [D, ones(length(D),1)*2], [E, ones(length(E),1)*3]);

GED_results_AB = griddata(xVals_AB,yVals_AB,GED_AB, ...
	dataSet_AB(:,1),dataSet_AB(:,2), 'nearest');

MED_results_AB = griddata(xVals_AB,yVals_AB,MED_AB, ...
	dataSet_AB(:,1),dataSet_AB(:,2), 'nearest');

MAP_results_AB = griddata(xVals_AB,yVals_AB,MAP_AB, ...
	dataSet_AB(:,1),dataSet_AB(:,2), 'nearest');

GED_results_CDE = griddata(xVals_CDE,yVals_CDE,GED_CDE, ...
	dataSet_CDE(:,1),dataSet_CDE(:,2), 'nearest');

MED_results_CDE = griddata(xVals_CDE,yVals_CDE,MED_CDE, ...
	dataSet_CDE(:,1),dataSet_CDE(:,2), 'nearest');

MAP_results_CDE = griddata(xVals_CDE,yVals_CDE,MAP_CDE, ...
	dataSet_CDE(:,1),dataSet_CDE(:,2), 'nearest');

confusionmat(dataSet_AB(:,3),GED_results_AB);
confusionmat(dataSet_AB(:,3),MED_results_AB);
confusionmat(dataSet_AB(:,3),MAP_results_AB);

cc=confusionmat(dataSet_CDE(:,3),GED_results_CDE);
dd=confusionmat(dataSet_CDE(:,3),MED_results_CDE);
ee=confusionmat(dataSet_CDE(:,3),MAP_results_CDE);


%% Test
A_test = gaussTransform(randn(N_C,2),mu_C,Sigma_A);
B_test = gaussTransform(randn(N_D,2),mu_D,Sigma_B);
C_test = gaussTransform(randn(N_C,2),mu_C,Sigma_C);
D_test = gaussTransform(randn(N_D,2),mu_D,Sigma_D);
E_test = gaussTransform(randn(N_E,2),mu_E,Sigma_E);

nn_dataSet_AB=vertcat([A_test, ones(length(A_test),1)*1], [B_test, ones(length(B_test),1)*2]);
nn_dataSet_CDE=vertcat([C_test, ones(length(C_test),1)*1], [D_test, ones(length(D_test),1)*2], [E_test, ones(length(E_test),1)*3]);

NN_results_AB = griddata(xVals_AB,yVals_AB,NN_class_AB, ...
	nn_dataSet_AB(:,1),nn_dataSet_AB(:,2), 'nearest');
NN5_results_AB = griddata(xVals_AB,yVals_AB,NN5_class_AB, ...
	nn_dataSet_AB(:,1),nn_dataSet_AB(:,2), 'nearest');

NN_results_CDE = griddata(xVals_CDE,yVals_CDE,NN_class_CDE, ...
	nn_dataSet_CDE(:,1),nn_dataSet_CDE(:,2), 'nearest');
NN5_results_CDE = griddata(xVals_CDE,yVals_CDE,NN5_class_CDE, ...
	nn_dataSet_CDE(:,1),nn_dataSet_CDE(:,2), 'nearest');

%% CONF?
confusionmat(nn_dataSet_AB(:,3),NN_results_AB)
confusionmat(nn_dataSet_AB(:,3),NN5_results_AB)
confusionmat(nn_dataSet_CDE(:,3),NN_results_CDE)
confusionmat(nn_dataSet_CDE(:,3),NN5_results_CDE)

