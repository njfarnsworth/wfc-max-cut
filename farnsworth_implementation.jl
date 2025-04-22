include("reusable_code.jl")
# using MaxCut

function main()

    file_name = "g1500-8861.dat"
    println("Graph Name: $file_name")
    max_cut = 0
    max_sets = Dict{Int, Int}()
    A = file_to_matrix(file_name)
    graph, edges = generate_graph(A) # graph generation process
    nodes = create_nodes_vector(graph, edges)
    println("# of Nodes: $(nv(graph))")
    println("# of Edges: $(ne(graph))")

    # viewgraph(graph)

    
    println("\nRunning 1000 iterations of WFC with parameters (200, 0.95, 200)\n")
    
    wfc_time = @elapsed begin
        for i in 1:1000
            cut, sets = wfc(nodes, edges, 200, .95, 200)
            if cut >= max_cut
                max_cut = cut
                max_sets = sets
            end
        end
    end
    println("WFC Time: $wfc_time")
    println("Max WFC Cut: $max_cut")
#=
   community = label_propagation_communities(graph)
   sol = genetics(community, max_sets, edges, 0.25, 100)

   println("This graph has $(length(community)) communities")
   println("Solution after genetic algorithm: ")
   for i in 1:length(sol)
    println(sol[i][2])
   end

     gw_time = @elapsed begin
        gw_cut, gw_set = maxcut(A)
    end

    println("\nGW Cut: $gw_cut")
    println("GW Time: $gw_time seconds")
    =#
end

main()