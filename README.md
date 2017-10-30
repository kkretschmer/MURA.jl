# MURA

[![Build Status](https://travis-ci.org/kkretschmer/MURA.jl.svg?branch=master)](https://travis-ci.org/kkretschmer/MURA.jl)
[![Coverage Status](https://coveralls.io/repos/github/kkretschmer/MURA.jl/badge.svg?branch=master)](https://coveralls.io/github/kkretschmer/MURA.jl?branch=master)

Generate [modified uniformly redundant arrays][mura] (MURAs) as
described in the publication below.  These arrays are used for [coded
aperture][ca] imaging.

[mura]: https://en.wikipedia.org/wiki/Modified_Uniformly_Redundant_Array
[ca]: https://en.wikipedia.org/wiki/Coded_aperture

> S. R. Gottesman and E. E. Fenimore. New family of binary arrays for
> coded aperture imaging. [Appl. Opt., 28:4344–4352, October
> 1989][ads]. doi:[10.1364/AO.28.004344](https://doi.org/10.1364/AO.28.004344).

[ads]: http://adsabs.harvard.edu/abs/1989ApOpt..28.4344G

Currently implements linear and square patterns and their decoding arrays.
Square patterns can be mosaicked.
