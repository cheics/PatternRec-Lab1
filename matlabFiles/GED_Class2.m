function [classNumber] = GED_Class2(dataPoint, ...
								meanA,covarA,n_A, ...
								meanB,covarB,n_B)
                            
[eigvec_A,eigval_A]=eig(covarA);
[eigvec_B,eigval_B]=eig(covarB);
whitetrans_A=inv(sqrt(eigval_A));
whitetrans_B=inv(sqrt(eigval_B));
                            
    function theDist = getDist(test_p, mu_test, whitetrans, eigvec)
        theDist=sqrt((test_p'-mu_test')'*eigvec*whitetrans*whitetrans*eigvec'...
            *(test_p'-mu_test'));
    end                        
    
    if getDist(dataPoint,meanB,whitetrans_B,eigvec_B)...
			>getDist(dataPoint,meanA,whitetrans_A,eigvec_A);
		classNumber=1;
	else
		classNumber=2;
	end
end