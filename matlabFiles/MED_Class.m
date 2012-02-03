function [classNumber] = MED_Class(dataPoint, varargin)
	classDistances=[];
	for k = 1:length(varargin) %Look at each mean
		currentMean=varargin{k};  
		squaredDist=sum((currentMean-dataPoint).^2);
		classDistances=[classDistances squaredDist];
	end
	[minimumDist, classNumber]=min(classDistances);
end