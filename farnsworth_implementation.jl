using Graphs, GraphRecipes, Plots, SimpleWeightedGraphs   # for plotting graphs
using StatsBase                     # for sample
using Combinatorics                 # for combinations
using Colors                        # to access RGB colors
using DataStructures                # for using stack and queue
using BenchmarkTools  

mutable struct Node
    key::Int
    degree::Int
    neighbors::Vector{Int}
    entropy::Float64
    set::Int
    edge_weights::Dict{Int, Float64} 
end


function make_random_graph(filename::AbstractString, num_nodes::Int, edge_multiplier::Real)
    filepath = joinpath(pwd(), "input_graphs", filename)
    open(filepath, "w") do io
        nodes = 1:num_nodes
        combos = collect(Combinatorics.combinations(nodes, 2))
        num_edges = trunc(Int, edge_multiplier * num_nodes)
        selected_edges = sample(combos, num_edges; replace=false)
        for (u, v) in selected_edges
            r = rand(1:10)
            println(io, "$u $v $r")
        end
    end
end

function read_data(filename)
    edges = Tuple{Int, Int, Float64}[]
    io = open(pwd() * "/input_graphs/" * filename, "r");
    for line in eachline(io)
        x, y, weight = [parse(Float64, ss) for ss in split(line)]
        t = Tuple([x, y, weight ])
        push!(edges, t)
    end
    return edges
end

function get_filename()
    i = 1
    while isfile(joinpath(pwd(), "input_graphs", "g$(i).txt"))
        i += 1
    end
    return "g$(i).txt"
end

function generate_file(num_nodes::Int = 7, edge_multiplier::Real = 1, filename::String = "nayda_graph.txt")
    if filename == ""
        filename = get_filename()
        make_random_graph(filename, num_nodes, edge_multiplier)
    end
    println("Using file: $filename")

    edges = read_data(filename)
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
        println("set: ", set)
        println("cut: ", cut)
        println("vertex: ", vertex, "\n\n")
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
    trials = 1 # declaring variables

    graph, edges = generate_file() # graph generation process
    nodes = create_nodes_vector(graph, edges)

    max_cut = 0

    for _ in 1:trials  
        cut, sets = wfc(graph, edges, 200, .95, 200, nodes)
        println(cut)
        println("\n\n\n\n", sets)
    end
end

main()