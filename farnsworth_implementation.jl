using Graphs, GraphRecipes, Plots, SimpleWeightedGraphs   # for plotting graphs
using StatsBase                     # for sample
using Combinatorics                 # for combinations
using Colors                        # to access RGB colors
using DataStructures                # for using stack and queue
using BenchmarkTools 
using MaxCut 
using Statistics

mutable struct Node
    key::Int
    degree::Int
    neighbors::Vector{Int}
    entropy::Float64
    set::Int
    edge_weights::Dict{Int, Float64} 
end

function matrix_generator()
    n = rand(3:10)
    A = zeros(n,n)
    for i in 1:size(A,1)
        for j in i:size(A,1)
            if i != j # diagonals must be 0
                A[i,j] = rand(1:20)
                A[j,i] = A[i,j] # for symmetry
            end
        end
    end
    return A
end

function matrix_to_edge(M)
    edge_list = []
    for i in 1:size(M,1)
        for j in i:size(M,1)
            if M[i, j] != 0
                push!(edge_list,(i,j,M[i,j]))
            end
        end
    end
    return(edge_list)
end

function generate_graph(M)
    edges = matrix_to_edge(M) # I wonder if I could just make the dictionary right here
    edge_list = Edge.(edges)  
    edge_dict = Dict{Tuple{Int,Int}, Float64}()
    for edge in edges
        edge_dict[(edge[1], edge[2])] = edge[3]
    end
    return SimpleGraph(edge_list), edge_dict # returns graph, edges
end

function create_nodes_vector(graph, edges)
    nodes = Vector{Node}()

    for key in 1:nv(graph)  # for each vertex in the graph
        deg = length(all_neighbors(graph, key))
        neighbor_weights = Dict{Int, Float64}()
        
        for neighbor in all_neighbors(graph, key)
            weight = get(edges, (key, neighbor), get(edges, (neighbor, key), 0.0))
            neighbor_weights[neighbor] = weight
        end
        
        new_node = Node(key, deg, all_neighbors(graph, key), 0, 0, neighbor_weights)
        push!(nodes, new_node)
    end
    
    return nodes

end

function propagate(v, nodes, sets, edges) 
    unpartitioned_neighbors = []
    propagatable_vertices = Vector{Node}()

    for node in nodes
        if node.key in v.neighbors && node.set == 0
            push!(unpartitioned_neighbors, node)
        end
    end

    for node in unpartitioned_neighbors
        can_propagate = true
        for neighbor_key in node.neighbors
            neighbor_node = get_node_by_key(nodes, neighbor_key)
            if neighbor_node === nothing || neighbor_node.set == 0
                can_propagate = false
                break
            end
        end
        if can_propagate
            push!(propagatable_vertices, node)
        end
    end

    for vertex in propagatable_vertices
        set, cut = propagate_cuts(vertex, nodes, edges)  # Use edge_dict instead of edges
        vertex.set = set
        push!(sets[set], vertex)
    end

    entropy_update(v, nodes, edges)
end

function collapse(v, nodes, sets, edges, temp, a)
    best_set, _ = propagate_cuts(v, nodes,edges)
    prob = rand()

    if (prob < exp(-a/temp))
        v.set = best_set == 1 ? 2 : 1
        push!(sets[v.set], v)
    else
        push!(sets[best_set], v)
        v.set = best_set
    end
end

function get_node_by_key(nodes::Vector{Node}, key::Int)
    for node in nodes
        if node.key == key
            return node
        end
    end
end

function calculate_cuts(edge_dict, sets)
    map = Dict{Int, Int}()
    for (index, set) in enumerate(sets)
        for v in set
            map[v.key] = index
        end
    end

    cut_edges = Set{Tuple{Int, Int, Float64}}()
    num_cuts = 0
    total_weight = 0.0

    for ((x, y), weight) in edge_dict
        if map[x] != map[y]
            num_cuts += 1
            push!(cut_edges, (x, y, weight))
            total_weight += weight
        end
    end

    return num_cuts, collect(cut_edges), total_weight
end

function entropy_update(v, nodes, edges)
    for key in v.neighbors
        neighbor = get_node_by_key(nodes, key)
        edge_weight = edges[(min(v.key, key), max(v.key, key))]
        neighbor.entropy += edge_weight
    end
end

function propagate_cuts(v, nodes, edges)
    max_1 = 0
    max_2 = 0

    for key in v.neighbors
        neighbor = get_node_by_key(nodes, key)
        edge_weight = edges[(min(v.key, key), max(v.key, key))]

        if neighbor.set == 1
            max_2 = max_2 + edge_weight
        end

        if neighbor.set == 2
            max_1 = max_1 + edge_weight
        end
    end
    return max_1 > max_2 ? (1, max_1) : (2, max_2)
end

function observe(nodes)
    unpartitioned_nodes = filter(node -> node.set ==0, nodes)
    max_entropy_node = isempty(unpartitioned_nodes) ? nothing : argmax(node -> node.entropy, unpartitioned_nodes)
    return(max_entropy_node)

end

function calculate_cuts(edge_dict, sets)
    map = Dict{Int, Int}()
    for (index, set) in enumerate(sets)
        for v in set
            map[v.key] = index
        end
    end

    cut_edges = Set{Tuple{Int, Int, Float64}}()
    num_cuts = 0
    total_weight = 0.0

    for ((x, y), weight) in edge_dict
        if map[x] != map[y]
            num_cuts += 1
            push!(cut_edges, (x, y, weight))
            total_weight += weight
        end
    end

    return total_weight
end


function wfc(graph, edges, temp, cooling, err_const, og_nodes)
    sets = [Set{Node}(), Set{Node}()] # two final partitioned sets

    unsorted_nodes = create_nodes_vector(graph, edges)
    nodes = sort(unsorted_nodes, by = x -> sum(values(x.edge_weights)), rev = true) # sorted nodes
    push!(sets[1], nodes[1]) # put the node with highest edge weights into the first set
    nodes[1].set = 1

    v = first(nodes) 
    x = propagate(v, nodes, sets, edges)

    while !isnothing(v)
        v = observe(nodes)
        if isnothing(v)
            break  
        end
        collapse(v, nodes, sets, edges, temp, err_const)
        propagate(v, nodes, sets, edges)
        temp *= cooling
    end

    return calculate_cuts(edges, sets), sets
end


function main()
    match = 0

    wfc_times = []
    gw_times = []

    for i in 1:1000
        A = matrix_generator()
        graph, edges = generate_graph(A) # graph generation process
        nodes = create_nodes_vector(graph, edges)
        max_cut = 0
        max_sets = []
        

        for _ in 1:1000
            wfc_time = @elapsed begin
                cut, sets = wfc(graph, edges, 100, .95, 100, nodes)
            end
            if cut >= max_cut
                max_cut = cut
                max_sets = sets
            end
            push!(wfc_times, wfc_time)
        end

        println("\n\n\nMax cut: ", max_cut)
        println("Max sets: ") 
        for set in max_sets
            print("Set: ")
            for node in set 
                print(node.key, " ")
            end
            println("")
        end

        gw_time = @elapsed begin
            GW_max_cut, max_partition = maxcut(A);
        end

        push!(gw_times, gw_time)


        println("GW Max Cut: ", GW_max_cut)

      
        if max_cut == GW_max_cut
            match += 1
        end
    end

    println("Cut matches: $match of 1000")
    println("Average WFC Time: $(mean(wfc_times))")
    println("Average GW Time: $(mean(gw_times))")
end

main()