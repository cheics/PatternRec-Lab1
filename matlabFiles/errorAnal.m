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

%% For nearest neighbor need to have training set...
dataSet_AB_shuf=dataSet_AB(randperm(length(dataSet_AB)), :);
dataSet_CDE_shuf=dataSet_CDE(randperm(length(dataSet_CDE)), :);

%% use half as training...
dataSet_AB_train=dataSet_AB_shuf(1:length(dataSet_AB)/2, :);
dataSet_AB_test=dataSet_AB_shuf(length(dataSet_AB)/2:end, :);
dataSet_CDE_train=dataSet_CDE_shuf(1:length(dataSet_CDE)/2, :);
dataSet_CDE_test=dataSet_CDE_shuf(length(dataSet_CDE)/2:end, :);

%% Establish classifiers
train_NN_class_AB=knnclassify(xyGrid_list_AB, ...
	dataSet_AB_train(:,1:2), ...
	dataSet_AB_train(:,3) ...
);
train_NN_class_AB=reshape(train_NN_class_AB,size(MED_AB,1),size(MED_AB,2));

train_NN_class_CDE=knnclassify(xyGrid_list_CDE, ...
	dataSet_CDE_train(:,1:2), ...
	dataSet_CDE_train(:,3) ...
);
train_NN_class_CDE=reshape(train_NN_class_CDE,size(MED_CDE,1),size(MED_CDE,2));

train_NN5_class_AB=knnclassify(xyGrid_list_AB, ...
	dataSet_AB_train(:,1:2), ...
	dataSet_AB_train(:,3), 5 ...
);
train_NN5_class_AB=reshape(train_NN5_class_AB,size(MED_AB,1),size(MED_AB,2));

train_NN5_class_CDE=knnclassify(xyGrid_list_CDE, ...
	dataSet_CDE_train(:,1:2), ...
	dataSet_CDE_train(:,3), 5 ...
);
train_NN5_class_CDE=reshape(train_NN5_class_CDE,size(MED_CDE,1),size(MED_CDE,2));

%% Test

NN_results_AB = griddata(xVals_AB,yVals_AB,train_NN_class_AB, ...
	dataSet_AB_test(:,1),dataSet_AB_test(:,2), 'nearest');
NN5_results_AB = griddata(xVals_AB,yVals_AB,train_NN5_class_AB, ...
	dataSet_AB_test(:,1),dataSet_AB_test(:,2), 'nearest');

NN_results_CDE = griddata(xVals_CDE,yVals_CDE,train_NN_class_CDE, ...
	dataSet_CDE_test(:,1),dataSet_CDE_test(:,2), 'nearest');
NN5_results_CDE = griddata(xVals_CDE,yVals_CDE,train_NN5_class_CDE, ...
	dataSet_CDE_test(:,1),dataSet_CDE_test(:,2), 'nearest');

%% CONF?
confusionmat(dataSet_AB_test(:,3),NN_results_AB);
confusionmat(dataSet_AB_test(:,3),NN5_results_AB);
confusionmat(dataSet_CDE_test(:,3),NN_results_CDE);
confusionmat(dataSet_CDE_test(:,3),NN5_results_CDE);

