namespace MARI.Optimization.NITRO

open MathNet.Numerics.LinearAlgebra.Double


module Derivatives=
        
        let DX = 1e-5

        type DIFFTYPE=
        |FORWARD
        |BACKWARD
        |CENTRAL
        
        let DEFAULT_DTYPE = CENTRAL
        
        let Diff(f:float->float,dx:float,dtype:DIFFTYPE)(x:float)=
            match dtype with
            |FORWARD -> (f(x+dx)-f(x))/dx
            |BACKWARD -> (f(x)-f(x-dx))/dx
            |CENTRAL -> (f(x+dx)-f(x-dx))/(2.0*dx)
            
        
        let Gradient(f:Vector -> float,dx:float,dtype:DIFFTYPE)(x:Vector):Vector=
            let g = Vector.Build.Dense(x.Count)  :?>Vector
            
            let f0 = f(x)
            match dtype with
            |FORWARD ->                        
                        for i in 0 .. x.Count-1 do
                            let y:Vector = x.Clone() :?>Vector
                            y.[i] <- y.[i]+dx
                            let a = f(y)
                            g.[i] <- (f(y)- f0)/dx
            |BACKWARD ->                                      
                        for i in 0 .. x.Count-1 do
                            let y = x.Clone()  :?>Vector
                            y.[i] <- y.[i]-dx
                            g.[i] <- (f0-f(y))/dx
            |CENTRAL ->                                                
                        for i in 0 .. x.Count-1 do
                            let yp = x.Clone()  :?>Vector
                            let ym = x.Clone()  :?>Vector
                            yp.[i] <- yp.[i]+dx
                            ym.[i] <- ym.[i]-dx
                            g.[i] <- (f(yp)-f(ym))/(2.0*dx)
            g
            
            
        let Jacobian(f:Vector->Vector,dx:float,dtype:DIFFTYPE)(x:Vector):Matrix=
            let f0 = f(x)
            let m=f0.Count
            let n = x.Count
            
            let J = Matrix.Build.Dense(n,m) :?> Matrix
            for j in 0 .. m-1 do
                let fj = fun z -> f(z).[j]
                J.SetColumn(j,Gradient(fj,dx,dtype)(x))
            
            J
            
        let Hessian(f:Vector -> float,dx:float,dtype:DIFFTYPE)(x:Vector):Matrix=
            let n = x.Count
            let h = Matrix.Build.Dense(n,n) :?>Matrix
                                                   
            for i in 0 .. n-1 do
                 let ei = Vector.Build.Dense(x.Count)  :?>Vector
                 ei.[i] <- dx
                                  
                 for j in i .. n-1 do       
                                                  
                    let ej = Vector.Build.Dense(x.Count)  :?> Vector
                    ej.[j] <- dx
             
                    let vx=
                        match dtype with
                        |FORWARD -> ( f((x+ei+ej)  :?>Vector )-f((x+ej)  :?>Vector) - (f((x+ei)  :?>Vector)-f(x)))/dx**2.0
                        |BACKWARD -> ( (f(x) - f((x-ei)  :?>Vector)) - (f((x-ej)  :?>Vector) - f((x-ei-ej)  :?>Vector)))/dx**2.0
                        |CENTRAL -> ((f((x+ei+ej) :?>Vector) - f((x-ei+ej) :?>Vector)) - (f((x+ei-ej)  :?>Vector) - f((x-ei-ej) :?>Vector)))/(4.0*dx**2.0) 
                                                   
                    h.[i,j] <- vx
            h


