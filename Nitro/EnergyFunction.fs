namespace MARI.Optimization.NITRO

open MathNet.Numerics.LinearAlgebra.Double

module EnergyFunction=
    open Derivatives    

    let E(f:Vector->float,h:Vector->Vector,lh:Vector,g:Vector -> Vector,lg:Vector,mu:float)
          (x:Vector,s:Vector):float=

          let gf = Gradient(f,DX,DEFAULT_DTYPE)(x)
          let Ah = Jacobian(h,DX,DEFAULT_DTYPE)(x)
          let Ag = Jacobian(g,DX,DEFAULT_DTYPE)(x)


          let a1 = (gf+Ah.Transpose()*lh+Ag.Transpose()*lg).InfinityNorm()
          let a2 = (Matrix.Build.DiagonalOfDiagonalArray(s.ToArray())*lg-mu).InfinityNorm()
          let a3 = h(x).InfinityNorm()
          let a4 = g(x).InfinityNorm()
          let a5 = (g(x)+s).InfinityNorm()

          [|a1;a2;a3;a4;a5|] |> Array.max

          
          
    let PHI(f:Vector->float,h:Vector->Vector,g:Vector->Vector,mu:float,nu:float)(x:Vector,s:Vector)=
        let barrier = s.ToArray() |> Array.map(fun x -> System.Math.Log(x)) |> Array.sum
        let normvec = Vector.Build.DenseOfArray(Array.concat([h(x).ToArray();(g(x)+s).ToArray()])).L2Norm()
        f(x) - mu * barrier + nu* normvec

