using MURA
using Primes
using Base.Test

const maxlength = 100
const randtests = 10

# Values from the publication
# 1989ApOpt..28.4344G, doi:10.1364_AO.28.004344
# Gottesman & Fenimore: New family of binary arrays for coded aperture imaging

# Table I. (p. 4)
@test linearlengths(74) == [5, 13, 17, 29, 37, 41, 53, 61, 73, 89]

# https://en.wikipedia.org/wiki/Quadratic_residue
@testset "quadratic residue" begin
    @test MURA.quadraticresidues(0) == []
    @test MURA.quadraticresidues(1) == [0]
    @test MURA.quadraticresidues(5) == [0, 1, 4]
    @test MURA.quadraticresidues(13) == [0, 1, 3, 4, 9, 10, 12]
    @test MURA.quadraticresidues(75) == [0, 1, 4, 6, 9, 16, 19, 21, 24,
                                         25, 31, 34, 36, 39, 46, 49, 51,
                                         54, 61, 64, 66, 69]
end

@testset "linear sequence" begin
    # Table I. (p. 4)
    #
    @testset "linear pattern" begin
        linearseq(L) =
            strip(replace(join(linearpattern(L)), r"([01]{5})", s"\1 "))
        patterns = Dict(
            05 => "01001",
            13 => "01011 00001 101",
            17 => "01101 00011 00010 11",
            29 => "01001 11101 00010 01000 10111 1001",
            37 => "01011 00101 11100 01000 01000 11110 10011 01",
            41 => "01101 10011 10000 01010 11010 10000 01110 01101 1",
            53 => (
                "01001 01101 11010 11100 00001 10011 00000 01110 10111 " *
                "01101 001"),
            61 => (
                "01011 10001 00111 11001 10100 10100 00001 01001 01100 " *
                "11111 00100 01110 1"),
            73 => (
                "01111 01011 00100 01011 00011 10100 00100 11110 01000 " *
                "01011 10001 10100 01001 10101 111"),
            89 => (
                "01101 10011 11000 01110 11100 10000 00101 01001 10101 " *
                "10101 10010 10100 00001 00111 01110 00011 11001 1011"),
        )
        for (L, P) in patterns
            @test linearseq(L) == P
        end
    end

    @testset "linear decoding" begin
        for L in linearlengths(maxlength)
            p, d = [L |> f |> symmshift for f in (linearpattern,
                                                  lineardecoding)]
            # apart from scaling, only one element differs
            @test (d .+ 1) .÷ 2 - p .|> abs |> sum == 1
            # test delta function like correlation
            @test sum(p .* d) == (L - 1) ÷ 2
            for i in 1:randtests
                s = rand(0:L - 1)
                @test sum(circshift(p, s) .* d) ==
                    (iszero(s) ? (L - 1) ÷ 2 : 0)
            end
        end
    end

    @testset "symmetry" begin
        A = "0011010101100"
        @test linearpattern(13) |> symmshift |> join == A
        for L in linearlengths(maxlength)
            for s in [L |> f |> symmshift for f in (linearpattern,
                                                    lineardecoding)]
                @test s == reverse(s)
            end
        end
    end

    @testset "argument checking" begin
        @test_throws ArgumentError linearpattern(6)
        @test_throws ArgumentError linearpattern(88)
    end
end # testset

@testset "square pattern" begin
    @testset "argument checking" begin
        @test_throws ArgumentError squarepattern(9)
        @test_throws ArgumentError linearpattern(57)
    end

    @testset "square mosaics" begin
        # The four n×n MURA aperture patterns from Fig. 5
        #
        patterns = Dict(
            5 => ["X..XXX..X",
                  ".XX.X.XX.",
                  ".XX.X.XX.",
                  "X..XXX..X",
                  ".........",
                  "X..XXX..X",
                  ".XX.X.XX.",
                  ".XX.X.XX.",
                  "X..XXX..X"],

            11 => ["X.XXX...X.XX.XXX...X.",
                   ".X...XXX.XX.X...XXX.X",
                   "X.XXX...X.XX.XXX...X.",
                   "X.XXX...X.XX.XXX...X.",
                   "X.XXX...X.XX.XXX...X.",
                   ".X...XXX.XX.X...XXX.X",
                   ".X...XXX.XX.X...XXX.X",
                   ".X...XXX.XX.X...XXX.X",
                   "X.XXX...X.XX.XXX...X.",
                   ".X...XXX.XX.X...XXX.X",
                   ".....................",
                   "X.XXX...X.XX.XXX...X.",
                   ".X...XXX.XX.X...XXX.X",
                   "X.XXX...X.XX.XXX...X.",
                   "X.XXX...X.XX.XXX...X.",
                   "X.XXX...X.XX.XXX...X.",
                   ".X...XXX.XX.X...XXX.X",
                   ".X...XXX.XX.X...XXX.X",
                   ".X...XXX.XX.X...XXX.X",
                   "X.XXX...X.XX.XXX...X.",
                   ".X...XXX.XX.X...XXX.X"]
        )

        a29 = "X..XXXX.X...X..X...X.XXXX..XXX..XXXX.X...X..X...X.XXXX..X"
        b29 = ".XX....X.XXX.XX.XXX.X....XX.X.XX....X.XXX.XX.XXX.X....XX."
        patterns[29] = [i == 'X' ? a29 : b29 for i in a29]
        patterns[29][29] = "........................................................."

        a59 = ".X...X.X.XX.XX...X....XX.....XXXXX..XXXX.XXX..X..X.X.XXX.XX.X...X.X.XX.XX...X....XX.....XXXXX..XXXX.XXX..X..X.X.XXX.X"
        b59 = "X.XXX.X.X..X..XXX.XXXX..XXXXX.....XX....X...XX.XX.X.X...X.XX.XXX.X.X..X..XXX.XXXX..XXXXX.....XX....X...XX.XX.X.X...X."
        patterns[59] = [i == 'X' ? a59 : b59 for i in a59]
        patterns[59][59] = "....................................................................................................................."

        # Go through the reference patterns and compare the computed mask
        # with the reference row-by-row. Use string comparison to easily
        # spot differences in test output
        for nelem in keys(patterns)
            mymask = squaremosaic(2nelem - 1)
            ref_mask = patterns[nelem]
            @test length(ref_mask[1]) in mosaiclengths(2nelem)
            for (n, ref_row) in enumerate(ref_mask)
                row = join([i == 1 ? "X" : "." for i in mymask[n, :]])
                @test row == ref_row
            end
        end
    end # testset

    @testset "square decoding" begin
        for L in primes(3, maxlength)
            p, d = [f(L) |> symmshift for f in (squarepattern,
                                                squaredecoding)]
            # apart from scaling, only one element differs
            @test (d .+ 1) .÷ 2 - p .|> abs |> sum == 1
            # test delta function like correlation
            @test sum(p .* d) == (L^2 - 1) ÷ 2
            for i in 1:randtests
                s = rand(0:L - 1, 2)
                @test sum(circshift(p, s) .* d) ==
                    (iszero(s) ? (L^2 - 1) ÷ 2 : 0)
            end
        end
    end

    @testset "symmetry" begin
        for L in primes(3, maxlength)
            for m in [f(L) |> symmshift for f in (squarepattern,
                                                  squaredecoding)]
                @test m == rot180(m)
            end
        end
    end
end # testset
