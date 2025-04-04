include("reusable_code.jl")

function main()

    max_cut = 0
    max_sets = []

    A = file_to_matrix("network.dat")
    graph, edges = generate_graph(A) # graph generation process
    
    wfc_time = @elapsed begin
        for i in 1:1000
            cut, sets = wfc(graph, edges, 200, .95, 200)
            if cut[3] >= max_cut
                max_cut = cut[3]
                max_sets = sets
            end
        end
    end
    println("WFC Time: $wfc_time")
    println("Max WFC Cut: $max_cut")
end

main()