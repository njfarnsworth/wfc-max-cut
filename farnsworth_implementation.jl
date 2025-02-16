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
    partitioned::Bool
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

function generate_file(num_nodes::Int = 7, edge_multiplier::Real = 1, filename::String = "g15.txt")
    if filename == ""
        filename = get_filename()
        make_random_graph(filename, num_nodes, edge_multiplier)
    end
    println("Using file: $filename")
    edge_list = Edge.(read_data(filename)) # 
    return SimpleGraph(edge_list), read_data(filename) # returns graph, edges
end

function create_nodes_vector(graph, edges)
    nodes = Vector{Node}()

    weights_map = Dict{Tuple{Int, Int}, Float64}()
    for edge in edges # add the edge weights to a dict 
        x, y, weight = edge
        weights_map[(x, y)] = weight
        weights_map[(y, x)] = weight  
    end
    
    for key in 1:nv(graph) # for each vertex in the graph, create the node struct
        deg = length(all_neighbors(graph, key))
        neighbor_weights = Dict{Float64, Float64}()
        for neighbor in all_neighbors(graph, key)
            neighbor_weights[neighbor] = weights_map[(key, neighbor)]
        end
        new_node = Node(key, deg, all_neighbors(graph, key), 0, 0, false, neighbor_weights)
        push!(nodes, new_node) 
    end
    return nodes
end

function calculate_cuts(edges, sets)
    map = Dict{Int, Int}() # this part just creates a map that determines which set each vertex belongs to
    for (index, set) in enumerate(sets)
        for v in set
            map[v.key] = index
        end
    end
    println(map)
    cut_edges = Set{Tuple{Int, Int, Float64}}()
    num_cuts = 0
    total_weight = 0
    for edge in edges
        x, y, weight = edge
        if map[x] != map[y]
            num_cuts += 1
            if !((y, x, weight) in cut_edges || (x, y, weight) in cut_edges)
                total_weight += weight
                push!(cut_edges, (x, y, weight))
            end
        end
    end
    return num_cuts, collect(cut_edges), total_weight
end


function propagate(v, nodes, sets, edges)
    unpartitioned_neighbors = []
    propagatable_vertices = []

    for node in nodes
        if node.key in v.neighbors && !node.partitioned # make these inter filters 
            push!(unpartitioned_neighbors, node)
        end
    end

    for node in unpartitioned_neighbors
        ct = 0
        for n in nodes
            if n.key in node.neighbors && !n.partitioned
                ct = ct + 1
            end
        end
        if ct == 0
            push!(propagatable_vertices, node)
        end
    end

    propagate_cuts(edges)
 ## breaks here, just trying to see which setting adding it to is better


end

function propagate_cuts(edges)
    edge_dict = Dict{Tuple{Int,Int}, Float64}()
    for edge in edges # want to do this forever eventually 
        edge_dict[(edge[1], edge[2])] = edge[3]
    end
    print(edge_dict)
end

function wfc(graph, edges, temp, cooling, err_const)
    sets = [Set{Node}(), Set{Node}()] # two final partitioned sets

    unsorted_nodes = create_nodes_vector(graph, edges)
    nodes = sort(unsorted_nodes, by = x -> sum(values(x.edge_weights)), rev = true) # sorted nodes
    push!(sets[1], nodes[1]) # put the node with highest edge weights into the first set
    nodes[1].partitioned = true
    nodes[1].set = 1

    v = popfirst!(nodes) # pop the node with highest edge weights
    x = propagate(v, nodes, sets, edges)

end


function main()
    trials = 1 # declaring variables

    graph, edges = generate_file() # graph generation process
    println("graph: ", graph)
    println("\n\n\nedges: ", edges)

    max_cut = 0

    for _ in 1:trials # do the wfc algorithm here 
        wfc(graph, edges, 200, .95, 200)
    end
    println("\n\n\n\n completed WFC")

end

main()