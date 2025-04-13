using Graphs, GraphRecipes, Plots, SimpleWeightedGraphs   
using StatsBase                   
using Combinatorics               
using Colors                       
using DataStructures               
using BenchmarkTools 
using Statistics
using Random

# NODE STRUCT 

mutable struct Node
    key::Int
    degree::Int
    neighbors::Vector{Int}
    entropy::Float64
    set::Int
    edge_weights::Dict{Int, Float64} 
end

# FILE READERS 

function tsplib_to_adj_matrix(filename::String) # turns certain file type into adjacency matrix 
    lines = readlines("instance_graphs/$filename")
    i0 = findfirst(x -> occursin("NODE_COORD_SECTION", x), lines)
    nodes = [parse.(Float64, tokens)[2:3] for l in lines[i0+1:end] 
             for tokens in [split(l)] if length(tokens) >= 3 && strip(l) != "EOF"]
    n = length(nodes)
    return [i == j ? 0 :
            floor(Int, sqrt((nodes[i][1]-nodes[j][1])^2 + (nodes[i][2]-nodes[j][2])^2) + 0.5)
            for i in 1:n, j in 1:n]
end


function tsp_adjacency_matrix(file_name, n) # turns different file type into adjacency matrix
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
 
function file_to_matrix(filename) # turns third file type into adjacency matrix .... need to consolidate all of these at some point 
    lines = readlines("/Users/naydafarnsworth/wfc-max-cut/community_detection_graphs/$filename")
    n = parse(Int, split(lines[1])[1])
    A = zeros(Float64, n, n)
        for line in lines[2:end]
            isempty(line) && continue
            parts = split(line)
            u = parse(Int, parts[1])
            v = parse(Int, parts[2])
            w = parse(Float64, parts[3])
            A[u, v] = w
            A[v,u] = w
        end
    return A
end


# GRAPH GENERATORS

function matrix_generator() # this generates complete graphs 
    n = rand(20:30)
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

# here --- make a random graph generator that generates a spanning tree and then makes it arbitraryily dense 

# OTHER

function matrix_to_edge(M) # makes a list of edges from an adjacency matrix 
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

function generate_graph(M) # creates a graph given an adjacency matrix 
    edges = matrix_to_edge(M) # I wonder if I could just make the dictionary right here
    edge_list = Edge.(edges)  
    edge_dict = Dict{Tuple{Int,Int}, Float64}()
    for edge in edges
        edge_dict[edge[1], edge[2]] = edge[3]
    end
    return SimpleGraph(edge_list), edge_dict # returns graph, edges
end

function create_nodes_vector(graph, edges) # creates a vector of nodes 
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

    return total_weight
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



function wfc(graph, edges, temp, cooling, err_const)
    final_sets = [Set{Node}(), Set{Node}()]
    final_cut = 0
   
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
            # propagate(v, nodes, sets, edges)
            temp *= cooling
        end
        cut = calculate_cuts(edges, sets)
        if cut > final_cut
            final_cut = cut
            final_sets = sets
        end

    return final_cut, final_sets
end


# Community detection functions 


function label_propagation_communities(G::Graph)
    coloring = _color_network(G)
    labeling = Dict(v => i for (i, v) in enumerate(vertices(G)))
    while !_labeling_complete(labeling, G)
        for (_color, nodes) in coloring
            for n in nodes
                _update_label(n, labeling, G)
            end
        end
    end
    clusters = Dict{Int, Set{Int}}()
    for (node, label) in labeling
        if !haskey(clusters, label)
            clusters[label] = Set{Int}()
        end
        push!(clusters[label], node)
    end
    return collect(values(clusters))
end

function greedy_color(G::Graph)
    colors = Dict{Int, Int}()
    for v in vertices(G)
        neighbor_colors = Set{Int}()
        for u in neighbors(G, v)
            if haskey(colors, u)
                push!(neighbor_colors, colors[u])
            end
        end
        color = 1
        while color in neighbor_colors
            color += 1
        end
        colors[v] = color
    end
    return colors
end

function _color_network(G::Graph)
    node_colors = greedy_color(G)
    coloring = Dict{Int, Set{Int}}()
    for (node, color) in node_colors
        if !haskey(coloring, color)
            coloring[color] = Set{Int}()
        end
        push!(coloring[color], node)
    end
    return coloring
end

function _update_label(n::Int, labeling::Dict{Int,Int}, G::Graph)
    label_counts = Dict{Int, Int}()
    for neighbor in neighbors(G, n)
        current_label = labeling[neighbor]
        label_counts[current_label] = get(label_counts, current_label, 0) + 1
    end
    if isempty(label_counts)
        return
    end
    max_label = first(keys(label_counts))
    max_count = label_counts[max_label]
    for (label, count) in label_counts
        if (count > max_count) || (count == max_count && label < max_label)
            max_label = label
            max_count = count
        end
    end
    labeling[n] = max_label
end

function _labeling_complete(labeling::Dict{Int,Int}, G::Graph)
    for n in vertices(G)
        label_counts = Dict{Int, Int}()
        for neighbor in neighbors(G, n)
            neighbor_label = labeling[neighbor]
            label_counts[neighbor_label] = get(label_counts, neighbor_label, 0) + 1
        end
        if isempty(label_counts)
            continue
        end
        max_label = first(keys(label_counts))
        max_count = label_counts[max_label]
        for (label, count) in label_counts
            if (count > max_count) || (count == max_count && label < max_label)
                max_label = label
                max_count = count
            end
        end
        if labeling[n] != max_label
            return false
        end
    end
    return true
end

function flip_community!(flipped_partition, community, sets)
    sets_in_comm = Set{Int}()
    for v in community
        push!(sets_in_comm, sets[v])
    end
    sets_vec = collect(sets_in_comm)
    mapping = Dict{Int, Int}()
    mapping[sets_vec[1]] = sets_vec[2]
    mapping[sets_vec[2]] = sets_vec[1]
    for v in community
        flipped_partition[v] = mapping[sets[v]]
    end
end

function selective_flips(communities, sets, edges, seq)
    flipped_partition = deepcopy(sets)
    for (i, community) in enumerate(communities)
        if seq[i] == 1
            flip_community!(flipped_partition, community, sets)
        end
    end
    new_cut_weight = 0.0
    for ((i, j), w) in edges
        if flipped_partition[i] != flipped_partition[j]
            new_cut_weight += w
        end
    end
    return flipped_partition, new_cut_weight, seq
end

function random_seq_flip(seq, p)
    for (i, val) in enumerate(seq)
        if rand() > p
            if val == 1
                seq[i] = 0
            else
                seq[i] = 1
            end
        end
    end
    return seq
end


function random_flips(communities, sets, edges, p)
    sequence = zeros(length(communities))
    flipped_partition = deepcopy(sets)
    for (i, community) in enumerate(communities)
        x = rand()
        if x >= p
            flip_community!(flipped_partition, community, sets)
            sequence[i] = 1
        end
    end
    new_cut_weight = 0.0
    for ((i, j), w) in edges
        if flipped_partition[i] != flipped_partition[j]
            new_cut_weight += w
        end
    end
    return flipped_partition, new_cut_weight, sequence
end


