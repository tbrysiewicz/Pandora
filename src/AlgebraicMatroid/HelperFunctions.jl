
export 
    rank_r_combination


function rank_r_combination(n :: Int, k :: Int, r ::Int)

    first :: Int64, sum :: Int64 = first_elt_comb(n,k,r)

    runningRank :: Int64 = r - sum

    comb :: Vector{Int64} = Vector{Int64}(undef, k)
    comb[1] = first


    for i in 1:(k - 1)

        first, sum = first_elt_comb(n - comb[i], k - i, runningRank)

        runningRank -= sum

        comb[i + 1] = first + comb[i]
        
    end

    return comb

end


function first_elt_comb(n :: Int, k :: Int, r ::Int)
    
    count :: Int = 1
    sum :: Int = binomial(n - 1, k - 1)

    while r > sum
        count += 1
        sum += binomial(n - count, k - 1)
    end
	
    return count, sum - binomial(n - count, k - 1)

end
