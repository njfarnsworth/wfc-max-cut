include("reusable_code.jl")
using MaxCut

function main()

    file_name = "g1000-7800.dat"

    println("Graph Name: $file_name")

    max_cut = 0
    max_sets = []

    A = file_to_matrix(file_name)
    graph, edges = generate_graph(A) # graph generation process

   #= println("# of Nodes: $(nv(graph))")
    println("# of Edges: $(ne(graph))")

    println("\nRunning 1000 iterations of WFC with parameters (200, 0.95, 200)\n")
    
    wfc_time = @elapsed begin
        for i in 1:1000
            cut, sets = wfc(graph, edges, 200, .95, 200)
            if cut >= max_cut
                max_cut = cut
                max_sets = sets
            end
        end
    end
    println("WFC Time: $wfc_time")
    println("Max WFC Cut: $max_cut")
=#
    gw_time = @elapsed begin
        gw_cut, gw_set = maxcut(A)
    end

    println("\nGW Cut: $gw_cut")
    println("GW Time: $gw_time seconds")

end

main()