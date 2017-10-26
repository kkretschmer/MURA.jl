module MURA

using Primes

export linearlengths, linearpattern, lineardecoding,
    squarepattern, squaremosaic, squaredecoding,
    mosaiclengths, symmshift

"""
    linearlengths(n)

Returns an array containing the valid MURA lengths up to the first one
greater or equal to `n`.  Valid lengths are primes of the class 4i +
1, where i is an integer (Gottesman & Fenimore 1989).
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

Returns an array containing the linear MURA sequence with length `L`
(Gottesman & Fenimore 1989, Eq. 11).
"""
function linearpattern(L::Integer)
    if !(L in linearlengths(L))
        throw(ArgumentError("L is not a valid length. See linearlengths()"))
    end
    q = quadraticresidues(L)
    return vcat([0], [i in q for i in 1:L - 1])
end

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

"""
    squarepattern(p)

Returns an array containing the square MURA pattern with `p` × `p`
elements (Gottesman & Fenimore 1989, Eq. 13 and 14).
"""
function squarepattern(p::Integer)
    if !isprime(p)
        throw(ArgumentError("p is not a valid length.  It must be prime."))
    end
    A = zeros(Int, p, p)
    q = quadraticresidues(p)
    A[2:end, 1] = 1
    for j in 2:p
        for i in 2:p
            if !((i - 1 in q) ⊻ (j - 1 in q))
                A[i, j] = 1
            end
        end
    end
    return A
end

"""
    squaremosaic(p)

Returns an array containing the square MURA mosaicked pattern with `p` ×
`p` elements (Gottesman & Fenimore 1989, Fig. 5).
"""
function squaremosaic(p::Integer)
    A = squarepattern(p ÷ 2 + 1)
    return hvcat((2, 2), A, A, A, A)[2:end, 2:end]
end

mosaiclengths(n) = primes(nextprime(ceil(Int, (n + 1) / 2))) * 2 - 1

"""
    squaredecoding(p)

Returns an array containing the decoding pattern for the square MURA
pattern with `p` × `p` elements (Gottesman & Fenimore 1989, Eq. 15).
"""
function squaredecoding(p::Integer)
    G = squarepattern(p) * 2 - 1
    G[1, 1] = 1
end

"""
    symmshift(A)

Rotate the array `A` to symmetrise a rectangular MURA pattern.
"""
symmshift(A) = circshift(A, size(A) .÷ 2)

end # module
