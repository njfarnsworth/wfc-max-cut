include("reusable_code.jl")
using MaxCut

function main()
    A = file_to_matrix("network.dat")
    gw_time = @elapsed begin
        GW_max_cut, max_partition = maxcut(A);
    end
    println("Cut: $GW_max_cut")
    println("Time: $gw_time")
end


main()