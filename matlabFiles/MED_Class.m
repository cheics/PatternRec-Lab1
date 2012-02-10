% varagin = variable # of arguments (in this cae means - mean a, mean b,
% etc.)
% classDistances is a vector of distances to the mean
% Calculates distance to each mean for a point
% Returns the minimum distance and class for the minimum

function [classNumber] = MED_Class(dataPoint, varargin)
	classDistances=[];
	for k = 1:length(varargin) %Look at each mean
		currentMean=varargin{k};  
		squaredDist=sum((currentMean-dataPoint).^2);
		classDistances=[classDistances squaredDist];
	end
	[minimumDist, classNumber]=min(classDistances);
end