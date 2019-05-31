function y=random_matrix1(M,N,s)

    y = zeros(M,N);
    
    for i=1:M
        ind = randperm(N);
        y(i,ind(1:s))=1; 
    end

end

