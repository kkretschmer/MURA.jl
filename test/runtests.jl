using MURA
using Base.Test

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

# Table I. (p. 4)
linearseq(L) = strip(replace(join(linearpattern(L)), r"([01]{5})", s"\1 "))
lineartest(L, seq) = @test linearseq(L) == seq

@testset "linear sequence" begin
    lineartest(05, "01001")
    lineartest(13, "01011 00001 101")
    lineartest(17, "01101 00011 00010 11")
    lineartest(29, "01001 11101 00010 01000 10111 1001")
    lineartest(37, "01011 00101 11100 01000 01000 11110 10011 01")
    lineartest(41, "01101 10011 10000 01010 11010 10000 01110 01101 1")
    lineartest(53,
               "01001 01101 11010 11100 00001 10011 00000 01110 10111 " *
               "01101 001")
    lineartest(61,
               "01011 10001 00111 11001 10100 10100 00001 01001 01100 " *
               "11111 00100 01110 1")
    lineartest(73,
               "01111 01011 00100 01011 00011 10100 00100 11110 01000 " *
               "01011 10001 10100 01001 10101 111")
    lineartest(89,
               "01101 10011 11000 01110 11100 10000 00101 01001 10101 " *
               "10101 10010 10100 00001 00111 01110 00011 11001 1011")
end
