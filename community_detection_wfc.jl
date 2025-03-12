using Graphs

# Define Vertex struct
mutable struct Vertex
    key::Int
    degree::Int
    neighbors::Vector{Int}
    entropy::Float64
    set::Int
    edge_weights::Dict{Int, Float64} 
end

function tsp_file_to_matrix(filename)
    lines = readlines(filename)
    n, _ = parse.(Int, split(lines[1]))
    A = zeros(Int, n, n)
    for line in lines[2:end]
        isempty(strip(line)) && continue
        u, v, w = parse.(Int, split(line))
        A[u, v] = w; A[v, u] = w
    end
    A
end

# Convert an adjacency matrix to a list of edges (with weights)
function matrix_to_edge(M)
    edge_list = []
    for i in 1:size(M,1)
        for j in i:size(M,1)
            if M[i, j] != 0
                push!(edge_list, (i, j, M[i,j]))
            end
        end
    end
    return edge_list
end

# Generate a SimpleGraph and an edge dictionary from the matrix
function generate_graph(M)
    edge_data = matrix_to_edge(M)
    # Build a list of unweighted edges for the graph.
    edge_list = [Edge(i, j) for (i, j, _) in edge_data]
    edge_dict = Dict{Tuple{Int,Int}, Float64}()
    for (i, j, w) in edge_data
        edge_dict[(i, j)] = w
    end
    n = size(M, 1)
    g = SimpleGraph(n)
    for e in edge_list
        add_edge!(g, src(e), dst(e))
    end
    return g, edge_dict
end

# --- Community Detection via Label Propagation ---

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

# --- WFC Procedure and Vertex Creation ---

function create_vertices_vector(graph, edge_dict)
    vertices_vec = Vector{Vertex}()
    for key in 1:nv(graph)
        deg = length(all_neighbors(graph, key))
        neighbor_weights = Dict{Int, Float64}()
        for neighbor in all_neighbors(graph, key)
            weight = get(edge_dict, (key, neighbor), get(edge_dict, (neighbor, key), 0.0))
            neighbor_weights[neighbor] = weight
        end
        new_vertex = Vertex(key, deg, all_neighbors(graph, key), 0.0, 0, neighbor_weights)
        push!(vertices_vec, new_vertex)
    end
    return vertices_vec
end

function get_vertex_by_key(vertices::Vector{Vertex}, key::Int)
    for v in vertices
        if v.key == key
            return v
        end
    end
    return nothing
end

function propagate(v, vertices, sets, edge_dict)
    unpartitioned_neighbors = []
    propagatable = Vector{Vertex}()
    for vertex in vertices
        if vertex.key in v.neighbors && vertex.set == 0
            push!(unpartitioned_neighbors, vertex)
        end
    end
    for vertex in unpartitioned_neighbors
        can_propagate = true
        for neighbor_key in vertex.neighbors
            neighbor = get_vertex_by_key(vertices, neighbor_key)
            if neighbor === nothing || neighbor.set == 0
                can_propagate = false
                break
            end
        end
        if can_propagate
            push!(propagatable, vertex)
        end
    end
    for vertex in propagatable
        set, _ = propagate_cuts(vertex, vertices, edge_dict)
        vertex.set = set
        push!(sets[set], vertex)
    end
    entropy_update(v, vertices, edge_dict)
end

function collapse(v, vertices, sets, edge_dict, temp, a)
    best_set, _ = propagate_cuts(v, vertices, edge_dict)
    prob = rand()
    if prob < exp(-a/temp)
        v.set = best_set == 1 ? 2 : 1
        push!(sets[v.set], v)
    else
        push!(sets[best_set], v)
        v.set = best_set
    end
end

function propagate_cuts(v, vertices, edge_dict)
    sum_set1 = 0.0
    sum_set2 = 0.0
    for key in v.neighbors
        neighbor = get_vertex_by_key(vertices, key)
        if neighbor === nothing
            continue
        end
        w = edge_dict[(min(v.key, key), max(v.key, key))]
        if neighbor.set == 1
            sum_set2 += w
        elseif neighbor.set == 2
            sum_set1 += w
        end
    end
    return sum_set1 > sum_set2 ? (1, sum_set1) : (2, sum_set2)
end

function entropy_update(v, vertices, edge_dict)
    for key in v.neighbors
        neighbor = get_vertex_by_key(vertices, key)
        if neighbor !== nothing
            w = edge_dict[(min(v.key, key), max(v.key, key))]
            neighbor.entropy += w
        end
    end
end

function observe(vertices)
    unpartitioned = filter(v -> v.set == 0, vertices)
    return isempty(unpartitioned) ? nothing : argmax(v -> v.entropy, unpartitioned)
end

function wfc(graph, edge_dict, temp, cooling, err_const, _)
    sets = [Set{Vertex}(), Set{Vertex}()]
    unsorted = create_vertices_vector(graph, edge_dict)
    vertices = sort(unsorted, by = v -> sum(values(v.edge_weights)), rev = true)
    push!(sets[1], vertices[1])
    vertices[1].set = 1
    v = first(vertices)
    propagate(v, vertices, sets, edge_dict)
    while true
        v = observe(vertices)
        if v === nothing
            break
        end
        collapse(v, vertices, sets, edge_dict, temp, err_const)
        propagate(v, vertices, sets, edge_dict)
        temp *= cooling
    end
    # Compute cut info (number of cut edges, list of cut edges, total weight)
    map_assignment = Dict{Int,Int}()
    for (i, s) in enumerate(sets)
        for v in s
            map_assignment[v.key] = i
        end
    end
    cut_edges = Set{Tuple{Int, Int, Float64}}()
    num_cuts = 0
    total_weight = 0.0
    for ((i, j), w) in edge_dict
        if map_assignment[i] != map_assignment[j]
            num_cuts += 1
            push!(cut_edges, (i, j, w))
            total_weight += w
        end
    end
    return (num_cuts, collect(cut_edges), total_weight), sets
end

# --- Wrap the Process as a Function ---
function run_max_cut(graph, full_edge_dict)

    # Detect communities.
    communities = label_propagation_communities(graph)
    
    # Global partition mapping: original vertex => set assignment.
    global_partition = Dict{Int, Int}()
    
    # Process each community.
    for i in 1:length(communities)
        community_vertices = sort(collect(communities[i]))
        # Create the induced subgraph.
        temp_result = induced_subgraph(graph, community_vertices)
        # Check if induced_subgraph returns a tuple; if so, extract the graph and mapping.
        if temp_result isa Tuple
            g_sub, mapping = temp_result
        else
            g_sub = temp_result
            mapping = community_vertices
        end

        # Build an edge dictionary for the subgraph.
        sub_edge_dict = Dict{Tuple{Int,Int}, Float64}()
        for e in edges(g_sub)
            sub_edge_dict[(src(e), dst(e))] = 1.0  # using a default weight; adjust as needed
        end
        
        # Run the WFC on the subgraph.
        (cut_info, sets) = wfc(g_sub, sub_edge_dict, 200, 0.95, 200, community_vertices)
        
        # Map the local partition back to the original graph.
        for s in sets
            for v in s
                orig_key = mapping[v.key]
                global_partition[orig_key] = v.set
            end
        end
    end
    
    # Compute the global cut weight using the full edge dictionary.
    global_cut_weight = 0.0
    for ((i, j), w) in full_edge_dict
        if global_partition[i] != global_partition[j]
            global_cut_weight += w
        end
    end
    return global_cut_weight, global_partition
end

function file_to_matrix(filename)
    lines = readlines("instance_graphs/$filename")
    n = parse(Int, split(lines[1])[1])
    A = zeros(Int, n, n)
        for line in lines[2:end]
            isempty(line) && continue
            u, v, w = parse.(Int, split(line))
            A[u, v] = w
            A[v,u] = w
        end
    return A
end


# --- Main Routine: Run 1000 Iterations and Choose the Best ---
function main()

    file_name = "mannino_k48.txt"
    A = file_to_matrix(file_name)
    g, e = generate_graph(A)
    best_cut = -Inf
    best_partition = Dict{Int, Int}()

    start_time = time()
    for i in 1:1000
        cut, partition = run_max_cut(g,e)
        if cut > best_cut
            best_cut = cut
            best_partition = partition
        end
    end

    elapsed = time() - start_time
    println("Max cut:", best_cut)
    println("Time: $elapsed")
end

main()
