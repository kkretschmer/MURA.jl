module MURA

using Primes

export linearlengths, linearpattern, lineardecoding

"""
    linearlengths(n)

Returns an array containing the valid MURA lengths up to the first one
greater or equal to `n`.
"""
function linearlengths(n)
    result = Int[]
    L = 5
    finished = false
    while !finished
        if isprime(L)
            append!(result, L)
            if L >= n
                finished = true
            end
        end
        L += 4
    end
    return result
end

"""
    quadraticresidues(n)

Returns an array containing the quadratic residues modulo `n`.
"""
quadraticresidues(n::Integer) = unique(sort((1:n).^2 .% n))

"""
    linearpattern(L)

Returns an array containing the linear MURA sequence with length `L`.
"""
    vcat([0], [i in quadraticresidues(L) for i in 1:L - 1])
function linearpattern(L::Integer)
"""
    lineardecoding(L)

Returns an array containing the decoding pattern for the linear MURA
sequence with length `L` (Gottesman & Fenimore 1989, Eq. 12).
"""
function lineardecoding(L::Integer)
    G = linearpattern(L) * 2 - 1
    G[1] = 1
    return G
end

end

end # module
