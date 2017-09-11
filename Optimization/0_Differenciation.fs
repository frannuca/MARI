namespace MARI.Optimization

open MathNet.Numerics.LinearAlgebra


module Differenciattion=
        
        type DIFFTYPE=
        |FORWARD
        |BACKWARD
        |CENTRAL
        
        
        let Diff(f:float->float,dx:float,dtype:DIFFTYPE)(x:float)=
            match dtype with
            |FORWARD -> (f(x+dx)-f(x))/dx
            |BACKWARD -> (f(x)-f(x-dx))/dx
            |CENTRAL -> (f(x+dx)-f(x-dx))/(2.0*dx)
            
        
        let Gradient(f:Vector<float> -> float,dx:float,dtype:DIFFTYPE)(x:Vector<float>)=
            let g = Vector<double>.Build.Dense(x.Count)
            
            let f0 = f(x)
            match dtype with
            |FORWARD ->                        
                        for i in 0 .. x.Count-1 do
                            let y = x.Clone()
                            y.[i] <- y.[i]+dx
                            let a = f(y)
                            g.[i] <- (f(y)- f0)/dx
            |BACKWARD ->                                      
                        for i in 0 .. x.Count-1 do
                            let y = x.Clone()
                            y.[i] <- y.[i]-dx
                            g.[i] <- (f0-f(y))/dx
            |CENTRAL ->                                                
                        for i in 0 .. x.Count-1 do
                            let yp = x.Clone()
                            let ym = x.Clone()
                            yp.[i] <- yp.[i]+dx
                            ym.[i] <- ym.[i]-dx
                            g.[i] <- (f(yp)-f(ym))/(2.0*dx)
            g
            
            
        let Jacobian(f:Vector<double>->Vector<double>,dx:float,dtype:DIFFTYPE)(x:Vector<float>)=
            let f0 = f(x)
            let m=f0.Count
            let n = x.Count
            
            let J = Matrix<double>.Build.Dense(n,m)
            for j in 0 .. m-1 do
                let fj = fun z -> f(z).[j]
                J.SetColumn(j,Gradient(fj,dx,dtype)(x))
            
            J
            
        let Hessian(f:Vector<float> -> float,dx:float,dtype:DIFFTYPE)(x:Vector<float>)=
            let n = x.Count
            let h = Matrix<double>.Build.Dense(n,n)
                                                   
            for i in 0 .. n-1 do
                 let ei = Vector<float>.Build.Dense(x.Count)
                 ei.[i] <- dx
                                  
                 for j in i .. n-1 do       
                                                  
                    let ej = Vector<float>.Build.Dense(x.Count)
                    ej.[j] <- dx
             
                    let vx=
                        match dtype with
                        |FORWARD -> ( f(x+ei+ej)-f(x+ej) - (f(x+ei)-f(x)))/dx**2.0
                        |BACKWARD -> ( (f(x) - f(x-ei)) - (f(x-ej) - f(x-ei-ej)))/dx**2.0
                        |CENTRAL -> ((f(x+ei+ej) - f(x-ei+ej)) - (f(x+ei-ej) - f(x-ei-ej)))/(4.0*dx**2.0) 
                                                   
                    h.[i,j] <- vx
            h