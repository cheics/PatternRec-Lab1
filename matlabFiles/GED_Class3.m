function [classNumber] = GED_Class3(dataPoint, ...
								meanC,covarC,n_C, ...
								meanD,covarD,n_D, ...
                                meanE,covarE,n_E)
                            
[eigvec_C,eigval_C]=eig(covarC);
[eigvec_D,eigval_D]=eig(covarD);
[eigvec_E,eigval_E]=eig(covarE);
whitetrans_C=inv(sqrt(eigval_C));
whitetrans_D=inv(sqrt(eigval_D));
whitetrans_E=inv(sqrt(eigval_E));
                            
    function theDist = getDist(test_p, mu_test, whitetrans, eigvec)
        theDist=sqrt((test_p'-mu_test')'*eigvec*whitetrans*whitetrans*eigvec'...
            *(test_p'-mu_test'));
    end

	comparePairs=[0 0 0];
	if getDist(dataPoint,meanD,whitetrans_D,eigvec_D)...
			>getDist(dataPoint,meanC,whitetrans_C,eigvec_C)
		comparePairs(1)=1;
	end
	if getDist(dataPoint,meanD,whitetrans_D,eigvec_D)...
			>getDist(dataPoint,meanE,whitetrans_E,eigvec_E)
		comparePairs(2)=1;
	end
	if getDist(dataPoint,meanE,whitetrans_E,eigvec_E)...
			>getDist(dataPoint,meanC,whitetrans_C,eigvec_C)
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