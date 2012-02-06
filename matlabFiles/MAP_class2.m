function [classNumber] = MAP_class2(dataPoint, ...
								meanA,covarA,n_A, ...
								meanB,covarB,n_B)

	function dist = getDist(dataPoint_t, mean_x,covar_x)
		dist= transpose(dataPoint_t'-mean_x')*inv(covar_x)*(dataPoint_t'-mean_x');
	end
	
	decision=getDist(dataPoint, meanB,covarB) ...
		-getDist(dataPoint, meanA,covarA) ...
		-2*log(n_B/n_A) ...
		-log(det(covarA)/det(covarB));
	
	if decision > 0
		classNumber=1;
	else
		classNumber=2;
	end


end