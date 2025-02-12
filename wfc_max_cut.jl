"""
This is the skeleton (that is, the starter code) for Pioneer Summer 2024 Mini-project 1.
"""

using Graphs, GraphRecipes, Plots, SimpleWeightedGraphs   # for plotting graphs
using StatsBase                     # for sample
using Combinatorics                 # for combinations
using Colors                        # to access RGB colors
using DataStructures                # for using stack and queue
using BenchmarkTools                # for assesing performance


mutable struct Node
    key::Int
    degree::Int
    neighbors::Vector{Int}
    entropy::Float64
    set::Int
    partitioned::Bool
    edge_weights::Dict{Int, Float64} 
end

"Write a new random graph to a file, where each line is a pair of node keys (e.g. 2 4)."
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

"Read a graph from a file, where each line is a pair of node keys (e.g. 2 4)"
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

"Display the given graph."
function viewgraph(graph, nodes)  
    p = graphplot(graph,
        names = 1:nv(graph),
        fontsize = 14,
        nodelabeldist = 5, 
        nodelabelangleoffset = π/4,
        markershape = :circle,
        markersize = 0.15,
        markerstrokewidth = 2,
        markerstrokecolor = :gray,
        edgecolor = :gray,
        linewidth = 2,
        curves = true
    )
    display(p)
end


function get_filename()
    i = 1
    while isfile(joinpath(pwd(), "input_graphs", "g$(i).txt"))
        i += 1
    end
    return "g$(i).txt"
end
    

"Determine the number of happy edges in the given graph (represented by its nodes). - for vertex coloring "
function count_happy_edges(nodes)
    nodecolors = [n.color for n in nodes]
    num_happies = 0
    for node in nodes
        for nbr in node.neighbors 
            num_happies += (node.color != nodecolors[nbr])
        end
    end
    return num_happies ÷ 2
end

"Creates an array of Node structs as a way of representing the given graph."
function create_nodes_vector(graph, edges)
    nodes = Vector{Node}()
    weights_map = Dict{Tuple{Int, Int}, Float64}()

    for edge in edges
        x, y, weight = edge
        weights_map[(x, y)] = weight
        weights_map[(y, x)] = weight  
    end
    

    for key in 1:nv(graph) # this loop is creating the nodes
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



function report_wfc(graph, edges, nodes)
    length(edges) < 63 && viewgraph(graph, nodes) 
end

############################################
pq = PriorityQueue{Node, Float64}()
sets = [Set{Node}(), Set{Node}()]
unpartitioned_vertices = []

function wfc_mc(graph, edges, temp, p, a)
    global sets
    global unpartitioned_vertices

    nodes = create_nodes_vector(graph, edges)
    println("IN WFC, NODES: ", nodes, "\n")
    unpartitioned_vertices = sort(nodes, by = x -> sum(values(x.edge_weights)), rev = true) # list of Nodes

    push!(sets[1], unpartitioned_vertices[1]) # puts the highest degree vertex into set 1
    unpartitioned_vertices[1].partitioned = true
    unpartitioned_vertices[1].set = 1
    vertex = popfirst!(unpartitioned_vertices)

    propagate(vertex, nodes) # good ish to here

    while length(unpartitioned_vertices) != 0
    
        vertex = observe()
        if isnothing(vertex)
            break
        end

        collapse(vertex, nodes, temp, a)
        propagate(vertex, nodes)
        temp *= p
    end

    return sets, nodes
    
end

function observe()
    global pq
    global unpartitioned_vertices

    vertex, entropy = dequeue_pair!(pq)
    if entropy == 0
        for v in unpartitioned_vertices
            push!(sets[1], v)
        end
        return nothing
    end


    return vertex
end

function collapse(vertex, nodes, temp, a)
    global unpartitioned_vertices
    global sets

    set1_count = 0
    set2_count = 0
    
    #finding number of neighbors in each set
    for v in vertex.neighbors
        if nodes[v].set == 1
            set1_count += nodes[v].edge_weights[vertex.key]
        end
    end

    set2_count = vertex.entropy - set1_count

    #finding set that resulted in more cuts
    if set1_count > set2_count
        good_set = 2
    elseif set1_count < set2_count
        good_set = 1
    else
        good_set = rand(1:2)
    end

    bad_set = 3 - good_set

    #Simulated Annealing
    prob = rand()

    if (prob < exp(-a/temp))
        push!(sets[bad_set], vertex)
        vertex.set = bad_set
    else
        push!(sets[good_set], vertex)
        vertex.set = good_set
    end


    vertex.partitioned = true
    filter!(x -> x != vertex, unpartitioned_vertices)

end

function propagate(vertex, nodes) # rewrite to partition neighbor nodes who have no unpartitioned neighbors
    global pq
    global sets

    up_neighbors = [n for n in vertex.neighbors if nodes[n].partitioned == false] #unpartitioned neighbors
    cont = false

    for n in up_neighbors
        set = nodes[nodes[n].neighbors[1]].set 
        for nn in nodes[n].neighbors
            if nodes[nn].partitioned == false # if any of the node's neighbors are unaartitioned, stop 
                cont = true
                break
            end
            if nodes[nn].set != set # if the neighbor is with the partitioned vertex, continue 
                cont = true
                break
            end
        end

        if cont
            continue
        end 

        push!(sets[3 - set], nodes[n]) # put the node in the opposite set  (kind of seems like this always happens)
        nodes[n].set = 3 - set
        nodes[n].partitioned = true
        filter!(x -> x != nodes[n], up_neighbors) # is he just removing one vertex or all but one ?
        filter!(x -> x != nodes[n], unpartitioned_vertices)
    end

    for n in up_neighbors
        nodes[n].entropy += nodes[n].edge_weights[vertex.key]
        if haskey(pq, nodes[n])
            delete!(pq, nodes[n])
        end
        enqueue!(pq, nodes[n], round(-(nodes[n].entropy), digits=0))
    end

    
    
end 
######################################################################

function calculate_cuts(edges, sets)

    map = Dict{Int, Int}() # this part just creates a map that determines which set each vertex belongs to
    for (index, set) in enumerate(sets)
        for v in set
            map[v.key] = index
        end
    end
  

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







#Benchmarking
function benchmark(graph, edges,  temp = 200, p = .95, a = 200)
    global pq = PriorityQueue{Node, Int}()
    global sets = [Set{Node}(), Set{Node}()]
    global unpartitioned_vertices = []
    subsets, edges = wfc_mc(graph, edges, temp, p, a)
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




#main function
function main(graph, edges, temp = 200, p = .95, a = 200, amt_trials = 10000) 
    # temperature changes, a fixed, p is what "cools" temp 
    max_cut = 0
    best_partition = []
    best_nodes = []

    for _ in 1:amt_trials
        sets, nodes = wfc_mc(graph, edges, temp, p, a)
        total_weight = calculate_cuts(edges, sets)
        if total_weight > max_cut
            max_cut = total_weight
            best_partition = deepcopy(sets)
            best_nodes = nodes
        end

        #reset global variables
        global pq = PriorityQueue{Node, Float64}()
        global sets = [Set{Node}(), Set{Node}()]
        global unpartitioned_vertices = []
    end

    println("---Wave Function Collapse 2---")
    println("Max Cut Weight: $max_cut")

end


graph, edges = generate_file()
println("Graph:", graph)
println("Edges:", edges)

main(graph, edges)
#@benchmark benchmark(graph, edges) 

        


