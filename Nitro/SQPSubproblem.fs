namespace MARI.Optimization.NITRO

open MathNet.Numerics.LinearAlgebra.Double


module SQPSubproblem=

    let Hxx(f:Vector->float)(x:Vector)=
        Derivatives.Hessian(f,Derivatives.DX,Derivatives.DEFAULT_DTYPE)


    let ComputeLagrange(Ah:Matrix,Ag:Matrix,S:Matrix,fgrad:Vector,mu:float)=
        let muvec = Vector.Build.Dense(Ah.ColumnCount+Ag.ColumnCount,mu)
        let Dv = Vector.Build.DenseOfArray(Array.concat([fgrad.ToArray();muvec.ToArray()]))

        let A = Matrix.Build.Dense(Ah.RowCount+S.RowCount,Ah.ColumnCount+Ag.ColumnCount)
        A.SetSubMatrix(0,0,Ah)
        A.SetSubMatrix(0,Ah.ColumnCount,Ag)
        A.SetSubMatrix(Ah.RowCount,Ah.ColumnCount,S)

        let r = ((A.Transpose()*A).Inverse()* A.Transpose()*Dv)
        r :?> Vector
    
    type PointDogleg={VCP:Vector;VN:Vector}


    let SearchDescendent(N:int,mint:float,maxt:float) (f:float -> bool)=
        let dt = (maxt-mint)/float(N)
        let mutable t = maxt;
            
        while t>=mint && not(f(t)) do
            t <- t - dt            
            
        if f(t) then
            Some(t)
        else 
            None


    let CauchyAndNewtonPoint(A:Matrix,h:Vector,g:Vector,s:Vector)=
        let V = Vector.Build.DenseOfArray(Array.concat([h.ToArray();(g+s).ToArray()]))
        let alpha_num = (A*V).L2Norm()**2.0
        let AAT = (A.Transpose()*A) :?>Matrix
        let alpha_den = V.ToRowMatrix()*(AAT*AAT)*V.ToColumnMatrix()
        let alpha = alpha_num/alpha_den.[0,0]
        
        let vCP = -alpha*A*V :?>Vector
        let vN = -1.0 * A*AAT.Inverse()*V :?>Vector
        {PointDogleg.VCP=vCP;PointDogleg.VN=vN}      

   

    let qfun(x:Vector,fgrad:Vector,mu:float,G:Matrix,Ns:int,Nx:int)=
        let mue = Vector.Build.Dense(Ns) + mu
        let V = Vector.Build.DenseOfArray(Array.concat [fgrad.ToArray();mue.ToArray()])

        let y = V.ToRowMatrix() * x + 0.5* x.ToRowMatrix()*G*x
        y.[0]
        
    let pred(d:Vector)=
        ()

    let NormalStep(A:Matrix,h:Vector,g:Vector,s:Vector,D:float,theta:float,tau:float)=
        
        let points = CauchyAndNewtonPoint(A,h,g,s)
        let vCP = points.VCP
        let vN = points.VN
        let N = 1000

        let solver = SearchDescendent(N,0.0,1.0)
        let t1 = solver(fun t -> 
                                   let sx = vN * t :?>Vector
                                   sx.L2Norm()<=tau*D &&
                                   vN.SubVector(g.Count,s.Count).L2Norm()>= -tau*0.5
                                   )

        let v =
            if t1.IsSome && t1.Value = 1.0 then
                vN
            else
                let t2= solver(fun t -> 
                                       let sx = (1.0-t)*vCP + vN * t :?>Vector
                                       sx.L2Norm()<=tau*D &&
                                       vN.SubVector(g.Count,s.Count).L2Norm()>= -tau*0.5
                                       )
                if t2.IsSome then
                    ((1.0-t2.Value)*vCP+t2.Value*vN ):?>Vector
                else
                    let t3 =  solver(fun t -> 
                                       let sx =  vCP * t :?>Vector
                                       sx.L2Norm()<=tau*D &&
                                       vN.SubVector(g.Count,s.Count).L2Norm()>= -tau*0.5
                                       )
                    (t3.Value*vCP) :?>Vector
      
        let vs = Matrix.Build.DiagonalOfDiagonalVector(s) * v.SubVector(g.Count,s.Count)
        v.SetSubVector(g.Count,s.Count,vs) 
        v
    let bisectionseach fw=
        let mutable tmax = 100.0
        let mutable tmin = 0.0
        let mutable t = (tmax+tmin)*0.5
                
        while System.Math.Abs(fw(t)-1.0)>1e-12 do 
            if fw(t)>1.0 then
                tmax <- t                        
            else 
                tmin <- t
                    
            t <- (tmax+tmin)*0.5
        
        t

    let TangentialStep(A:Matrix,Hxx:Matrix,Hss:Matrix,gradf:Vector,mu:float,vnormal:Vector,D:float)=
        let G = Matrix.Build.Dense(Hxx.RowCount+Hss.RowCount,Hxx.ColumnCount+Hss.ColumnCount)
        G.SetSubMatrix(0,0,Hxx)
        G.SetSubMatrix(Hxx.RowCount,Hxx.ColumnCount,Hss)
        let N = Hxx.RowCount+Hss.RowCount
        let P = Matrix.Build.DenseIdentity(N) - A*(A.Transpose()*A).Inverse() *A.Transpose()

        let muv = mu * Vector.Build.Dense(Hss.RowCount,1.0)
        let mutable w = Vector.Build.Dense(Hxx.RowCount+Hss.RowCount)
        let mutable r = Vector.Build.DenseOfArray(Array.concat([gradf.ToArray();muv.ToArray()])) + G*vnormal
        let mutable g = P*r
        let mutable p = -1.0 * g
        let aa = g.ToRowMatrix()*r
        let tol = 0.01*sqrt(aa.[0])


        
        let itermax = Hxx.RowCount-(A.ColumnCount-Hss.ColumnCount)
        let counter = ref 0
        while !counter< itermax do
            let c = p.ToRowMatrix()*G*p
            let fw  (t:float) = (w + t*p).L2Norm()/D

            if c.[0] <= 0.0 then
                                
                let mutable tmax = 100.0
                let mutable tmin = 0.0
                let mutable t = (tmax+tmin)*0.5
                
                let t=bisectionseach(fw)                
                w <- w+t*p   
                counter := itermax
            else
                let alpha = r.ToRowMatrix()*g/c.[0]
                w <- w + alpha.[0]*p
                if w.L2Norm() > D then
                    let t=bisectionseach(fw)                
                    w <- w+t*p    
                    counter := itermax
                else
                    r <- r + alpha.[0]*G*p
                    g <- P*r
                    let error = g.ToRowMatrix()*r.ToColumnMatrix()
                    if error.[0,0] < tol then
                        counter := itermax
                    
                    let aux = r.ToRowMatrix()*g
                    let beta = r.ToRowMatrix()*g/aux.[0]
                    p <- -g + beta.[0]*p
                    
                        
        w                
                

    let SQP(mu:float,eu:float,fgrad:Vector,hxx:Matrix,hss:Matrix,hk:Vector,gk:Vector,Ag:Matrix,Ah:Matrix,D:float,theta:float,tau:float)(xk:Vector,sk:Vector)=
        
        let S  = Matrix.Build.DenseOfDiagonalArray(sk.ToArray()) :?> Matrix
        let L = ComputeLagrange(Ah,Ag,S,fgrad,mu)
        let A = Matrix.Build.Dense(Ah.RowCount+S.RowCount,Ah.ColumnCount+Ag.ColumnCount)
        A.SetSubMatrix(0,0,Ah)
        A.SetSubMatrix(0,Ah.ColumnCount,Ag)
        A.SetSubMatrix(Ah.RowCount,Ah.ColumnCount,S)

        
        let Lg = L.SubVector(Ah.RowCount,Ag.RowCount) :?> SparseVector

        let Hss = Matrix.Build.Dense(sk.Count,sk.Count)

        for i in 0 .. Lg.Count do
            if Lg.[i]/sk.[i] > 0.0 then
                Hss.[i,i]<-Lg.[i]/sk.[i]
            else
                Hss.[i,i]<-mu/(sk.[i])**2.0
        

        let v=NormalStep(A :?> Matrix ,hk,gk,sk,D,theta,tau)
        let w = TangentialStep(A:?> Matrix,hxx,hxx,fgrad,mu,v,D)
        
        let dk = v + w

        //TODO: 
        
        failwith "Not implemented"

