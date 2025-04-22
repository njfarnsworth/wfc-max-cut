include("reusable_code.jl")

function main()
    file_name = "g2000-15402.dat"
    println("Graph Name: $file_name")
    A = file_to_matrix(file_name)
    graph, edges = generate_graph(A) 
    all_edges, internal_edges = edges_between_communities(graph)
    println("There are $all_edges edges in the graphs and $(all_edges-internal_edges) edges connecting communities")
end

main()