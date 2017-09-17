namespace MARI.Optimization.NITRO

open MathNet.Numerics.LinearAlgebra.Double


module NitroAlgorithm=
    let Optimize(f:Vector->float,h:Vector->Vector,g:Vector->Vector,eu:float,mu:float,etol:float)(x:Vector)=

        let aux_g = g(x)
        let aux_h = h(x)
        let lh = ref (Vector.Build.Dense(aux_h.Count,0.0))
        let lg = ref (Vector.Build.Dense(aux_g.Count,0.0))
        let tol= EnergyFunction.E(f,h,!lh :?>Vector, g,!lg :?>Vector ,0.0)
        let s = ref (Vector.Build.Dense(aux_g.Count,0.0))
        while tol(x,!s :?>Vector) > etol do
            ()
        ()

