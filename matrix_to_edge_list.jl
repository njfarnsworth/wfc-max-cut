A = [0 2 3.0; 2 0 1.5; 3.0 1.5 0]

function matrix_to_edge(M)
    edge_list = []
    for i in 1:size(M,1)
        for j in i:size(M,1)
            if A[i, j] != 0
                push!(edge_list,(i,j,A[i,j]))
            end
        end
    end
    return(edge_list)
end

println(matrix_to_edge(A))