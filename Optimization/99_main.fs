namespace MARI.Optimization

module Program =
    
    open MathNet.Numerics.LinearAlgebra
    
    open Differenciattion
    
    let f (x:Vector<float>) = (x.[0]-2.0)**2.0 + x.[1]**2.0
    
    let feq (x:Vector<float>) = Vector<float>.Build.DenseOfArray( [| 2.0*x.[0] - 10.0*x.[1];
                                                                     System.Math.Sin(x.[0]*3.0)|])
        
    
    [<EntryPoint>]
    let main args =
        let x0 = Vector.Build.DenseOfArray([|0.0;0.0|])
        
        let printmatrix dtype=
            let g = Gradient(f,1e-5,dtype)(x0)
            let h = Hessian(f,1e-5,dtype)(x0)
            let j = Jacobian(feq,1e-5,dtype)(x0)
            System.Console.WriteLine(g.ToString())
            System.Console.WriteLine(h.ToString())
            System.Console.WriteLine(j.ToString())
        
        printmatrix DIFFTYPE.CENTRAL
        printmatrix DIFFTYPE.FORWARD
        printmatrix DIFFTYPE.CENTRAL
                
                    
                    
        
        0