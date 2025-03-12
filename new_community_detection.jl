using Graphs



function create_nodes(graph)
    nodes = Vector{Node}()
    for key in 1:nv(graph) 
        new_node = Node(key, all_neighbors(graph, key), key)
        push!(nodes, new_node)
    end
    return nodes
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

function get_node_by_key(nodes::Vector{Node}, key::Int)
    for node in nodes
        if node.key == key
            return node
        end
    end
end
