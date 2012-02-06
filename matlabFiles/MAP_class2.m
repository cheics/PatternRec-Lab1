function [classNumber] = MAP_class2(dataPoint, ...
								meanA,covarA,n_A, ...
								meanB,covarB,n_B)

	function theDist = getDist(test_p, mu_test, covar_test)
		theDist=(transpose(test_p'-mu_test')) ...
					*(inv(covar_test)) ...
					*(test_p'-mu_test');
	end

	if getDist(dataPoint,meanB,covarB)-getDist(dataPoint,meanA,covarA) ...
			>2*log(n_B/n_A) + log(det(covarA)/det(covarB))
		classNumber=1;
	else
		classNumber=2;
	end
end