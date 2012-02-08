dataSet_AB=vertcat([A, ones(length(A),1)*1], [B, ones(length(B),1)*2]);
dataSet_CDE=vertcat([C, ones(length(C),1)*1], [D, ones(length(D),1)*2], [E, ones(length(E),1)*3]);

GED_results_AB = interp2(xVals_AB,yVals_AB,GED_AB, ...
	dataSet_AB(:,1),dataSet_AB(:,2), 'nearest');

MED_results_AB = interp2(xVals_AB,yVals_AB,MED_AB, ...
	dataSet_AB(:,1),dataSet_AB(:,2), '*nearest*');

MAP_results_AB = interp2(xVals_AB,yVals_AB,MAP_AB, ...
	dataSet_AB(:,1),dataSet_AB(:,2), '*nearest*');

GED_results_CDE = interp2(xVals_CDE,yVals_CDE,GED_CDE, ...
	dataSet_CDE(:,1),dataSet_CDE(:,2), '*nearest*');

MED_results_CDE = interp2(xVals_CDE,yVals_CDE,MED_CDE, ...
	dataSet_CDE(:,1),dataSet_CDE(:,2), '*nearest*');

MAP_results_CDE = interp2(xVals_CDE,yVals_CDE,MAP_CDE, ...
	dataSet_CDE(:,1),dataSet_CDE(:,2), '*nearest*');