function [classNumber] = MAP_class3(dataPoint, ...
								meanA,covarA,n_A, ...
								meanB,covarB,n_B, ...
								meanC,covarC,n_C)

	function theDist = getDist(test_p, mu_test, covar_test)
		theDist=(transpose(test_p'-mu_test')) ...
					*(inv(covar_test)) ...
					*(test_p'-mu_test');
	end

	comparePairs=[0 0 0];
	if getDist(dataPoint,meanB,covarB)-getDist(dataPoint,meanA,covarA) ...
			>2*log(n_B/n_A) + log(det(covarA)/det(covarB))
		comparePairs(1)=1;
	end
	if getDist(dataPoint,meanB,covarB)-getDist(dataPoint,meanC,covarC) ...
			>2*log(n_B/n_C) + log(det(covarC)/det(covarB))
		comparePairs(2)=1;
	end
	if getDist(dataPoint,meanC,covarC)-getDist(dataPoint,meanA,covarA) ...
			>2*log(n_C/n_A) + log(det(covarA)/det(covarC))
		comparePairs(3)=1;
	end
	
	if comparePairs(1)==1 && comparePairs(3)==1
		classNumber=1;
	elseif comparePairs(1)==0 && comparePairs(2)==0
		classNumber=2;
	elseif comparePairs(2)==1 && comparePairs(3)==0
		classNumber=3;
	else
		classNumber=-1
	end
end