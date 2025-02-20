using MaxCut

W = [0 8 2 1 0; 
     8 0 3 2 0; 
     2 3 0 0 0; 
     1 2 0 0 1; 
     0 0 0 1 0];


function rewrite(W)
    open("/input_graphs/edges.txt", "w") do file
        for i in 1:size(W, 1)
            for j in i+1:size(W, 2)  # Avoid duplicate edges in an undirected graph
                if W[i, j] != 0
                    println(file, "$i, $j, $(W[i, j])")
                end
            end
        end
    end
println("Edges saved to edges.txt")
end



max_cut, max_partition = maxcut(W);






	
@show max_cut;
## max_cut = 14

@show max_partition;
## max_partition = ([1, 3, 4], [2, 5])