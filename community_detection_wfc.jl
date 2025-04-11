include("reusable_code.jl")




function main()

    file_name = "network.dat"
    A = file_to_matrix(file_name)
    global  g, e = generate_graph(A)
    best_cut = -Inf
    global best_partition = Dict{Int, Int}()

    cd_time = @elapsed begin
        global communities = label_propagation_communities(g)
        println("\nThis graph has $(length(communities)) communities.")
        global_cut_weight = 0.0

        println("\nRunning 1000 iterations of WFC + CD with parameters (200, 0.95, 200)\n")

        for _ in 1:1000
            global_partition = Dict{Int, Int}()
            for i in 1:length(communities)
                community_vertices = sort(collect(communities[i]))
                temp_result = induced_subgraph(g, community_vertices)
                if temp_result isa Tuple
                    g_sub, mapping = temp_result
                else
                    g_sub = temp_result
                    mapping = community_vertices
                end
                sub_edge_dict = Dict{Tuple{Int,Int}, Float64}()
                for edge in edges(g_sub)
                    u_sub = src(edge)
                    v_sub = dst(edge)
                    u_orig = mapping[u_sub]
                    v_orig = mapping[v_sub]
                    key = (min(u_orig, v_orig), max(u_orig, v_orig))
                    weight = e[key]
                    sub_edge_dict[(u_sub, v_sub)] = weight
                end
                _, sets = wfc(g_sub, sub_edge_dict, 200, .95, 200)  
                for s in sets
                    for v in s
                        orig_key = mapping[v.key]
                        global_partition[orig_key] = v.set
                    end
                end
            end
            global_cut_weight = 0.0
            for ((i, j), w) in e
                if global_partition[i] != global_partition[j]
                    global_cut_weight += w
                end
            end
            if global_cut_weight > best_cut
                best_cut = global_cut_weight
                best_partition = global_partition
            end
        end
    end

    println("WFC + CD Time: $cd_time seconds")
    println("WFC + Community detection Max Cut: ", best_cut)
    

end

main()
