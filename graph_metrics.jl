using Graphs, SparseArrays, LinearAlgebra

function tsp_adjacency_matrix(file_name, n)
    data = read("instance_graphs/$file_name", String)
    weights = parse.(Int, split(data))

    matrix = zeros(Int, n, n)
    index = 1
    for i in 1:n
        for j in 1:i  
            matrix[i, j] = weights[index]
            matrix[j, i] = weights[index] 
            index += 1
        end
    end
    return matrix
end



function compute_graph_invariants(g::AbstractGraph)
    println("=== Graph Invariants ===")
    
    # Basic properties
    println("Number of Nodes (Order): ", nv(g))
    println("Number of Edges (Size): ", ne(g))
    println("Graph Density: ", 2 * ne(g) / (nv(g) * (nv(g) - 1)))

    # Degree-based invariants
    degrees = degree(g)
    println("Maximum Degree: ", maximum(degrees))
    println("Minimum Degree: ", minimum(degrees))
    println("Average Degree: ", sum(degrees) / nv(g))
    
    # Connectivity
    println("Is Connected: ", is_connected(g))
    println("Diameter (if connected): ", is_connected(g) ? diameter(g) : "Graph is disconnected")
    
    # Spectral properties
    laplacian_matrix = sparse(Diagonal(degrees) - adjacency_matrix(g))
    laplacian_eigenvalues = eigvals(Matrix(laplacian_matrix))
    
   # println("Laplacian Eigenvalues: ", laplacian_eigenvalues)
    println("Spectral Gap (λ2 - λ1): ", laplacian_eigenvalues[2] - laplacian_eigenvalues[1])
    
    # Clustering coefficient (approximation)
    triangle_count = sum(neighbors(g, v) ∩ neighbors(g, w) != [] for v in vertices(g) for w in neighbors(g, v)) ÷ 2
    println("Number of Triangles: ", triangle_count)
    

end


graph = Graph(tsp_adjacency_matrix("pw01_100_0.txt", 100))
compute_graph_invariants(graph)
