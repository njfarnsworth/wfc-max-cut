include("reusable_code.jl")
include("community_detection_wfc.jl")


function genetics(communities, best_partition, e, p, gen)
    sols = []
    for _ in 1:1000
        sets, cut, sequence = random_flips(communities, best_partition, e, p)
        push!(sols, (sets, cut, sequence))
    end
    sorted_entries = sort(sols, by = x -> x[2], rev = true)
    current_top_solutions = sorted_entries[1:100]
    for _ in 1:gen
        children = []
        for _ in 1:100
            j = rand(1:length(communities)-1)
            _, _, seq1 = rand(current_top_solutions)
            _, _, seq2 = rand(current_top_solutions)
            new_seq = vcat(seq1[1:j], seq2[j+1:end])
            child_sequence = random_seq_flip(new_seq, 0.02)
            child = selective_flips(communities, best_partition, e, child_sequence)
            push!(children, child)
        end
        sorted_children = sort(children, by = x -> x[2], rev = true)
        current_top_solutions = sorted_children[1:min(100, length(sorted_children))]
    end

    return current_top_solutions[1:10]
end


function main()

    println("\n Genetic Algorithm Parameters:")
    println("\t * 1000 parents in first generation")
    println("\t * 25% chance of community flip for first parent generation")
    println("\t * Use top 100 solutions to generate 1000 children, keep 100, repeat")
    println("\t * 2% mutation rate for children")
    println("\t * 100 generations total\n")

    genetics_time = @elapsed begin 
        sols = genetics(communities, best_partition, e, 0.25, 100)
    end

    top_sol = sols[1][2]

    println("Additional time for genetic algorithm: $genetics_time")
    println("Top solution after genetics algorithm: $top_sol")

    
    # we have a solution - let 0 represent an original community and 1 represent a flipped community. Create a list of 0s representing the original communities
    
    # for the first generation of solutions, we create (say 100) new community possibilities by performing random flips. 
    # For our sake, let's have a probability and then go trhough and apply it to each community. Then, take the cut of all these new ones and keep the top 10% of the results

    # now, we enter a loop in which we breed "better answers", take theirr "children" and "breed them"
    # generate random couples of 
end

main()