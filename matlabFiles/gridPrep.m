function [xVals, yVals, classGrid] = gridPrep(gridSize, varargin)
	% Input gridsize and the (2D) data sets
	minX_vals=[];
	minY_vals=[];
	maxX_vals=[];
	maxY_vals=[];
	% For each dataset find the bounds
	for k = 1:length(varargin)
		xVals=varargin{k}(:, 1);           % Cell array indexing
		yVals=varargin{k}(:, 2);
		minX_vals = [minX_vals min(xVals)];
		minY_vals = [minY_vals min(yVals)];
		maxX_vals = [maxX_vals max(xVals)];
		maxY_vals = [maxY_vals max(yVals)];
	end
	
	xVals = min(minX_vals)-1:gridSize:max(maxX_vals)+1;
	yVals = min(minY_vals)-1:gridSize:max(maxY_vals)+1;
	classGrid = zeros(length(yVals),length(xVals));

end
