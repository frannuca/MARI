namespace MARI.Optimization


module InteriorPoint=
    open MathNet.Numerics.LinearAlgebra
    
    
    module NITRO =
        let E (xk:Vector<float>, sk:Vector<float>, mu:float)=
                1.0
     
        let SearchDescendent(N:int,mint:float,maxt:float) (f:float -> bool)=
            let dt = (maxt-mint)/float(N)
            let mutable t = maxt;
            
            while t>=mint && not(f(t)) do
                t <- t - dt
            
            if f(t) then
                Some(t)
            else 
                None
                
            
        let Dogleg (mf:Matrix<float> -> float,Ah:Matrix<float>,Ag:Matrix<float>, h:Vector<float>,g:Vector<float>,x:Vector<float>,s:Vector<float>,
                    theta:float,zeta:float,Delta:float,tau:float) =
                    
            let S = Matrix<float>.Build.DenseOfDiagonalVector(s)
            let Ak = Matrix<float>.Build.Dense(Ah.RowCount+S.RowCount,Ah.ColumnCount+Ag.ColumnCount)
            Ak.SetSubMatrix(0,0,Ah)
            Ak.SetSubMatrix(0,Ah.ColumnCount,Ag)
            Ak.SetSubMatrix(Ah.RowCount,Ah.ColumnCount,S)
            
            let AA = (Ak.Transpose()*Ak)
            let AA2 = AA * AA
            
            let hgs = Vector<float>.Build.DenseOfArray(Array.concat([h.AsArray();(g+s).AsArray()]))
            let alpha_num =  (Ak * hgs)
            let alpha_den = hgs.ToRowMatrix()*AA2*hgs.ToColumnMatrix()
            let alpha = alpha_num.L2Norm()/(alpha_den.[0,0])
            
            //Cauchy point
            let vCP = -alpha * Ak * hgs.ToColumnMatrix()
            
            //Normal direction:
            let Vn = -Ak * AA.Inverse() * hgs.ToColumnMatrix()
            
            let normVn = Vn.L2Norm()
            
            let th1 = min 1.0 (zeta*Delta/normVn)
            
            
            
            let isSvalid (vx:Matrix<float>) = vx.L2Norm() <= zeta*Delta && 
                                              vx.SubMatrix(x.Count,s.Count,0,1).AsColumnMajorArray() 
                                              |> Array.exists(fun xs -> xs >= -tau*0.5)
            
            let ft1 =   (fun  t -> (1.0-t)*vCP+t*Vn)>>isSvalid
            let t1 = SearchDescendent(1000,0.0,1.0) ft1
                       
            
            
            
            let v = 
                match t1 with
                |Some(t1) when t1=1.0 ->  Vn
                |Some(t1) ->
                            let ft2 = (fun t -> (1.0-t)*vCP+t*Vn) >> isSvalid
                            let t2 = SearchDescendent(1000,0.0,1.0)(ft2)
                            let vDL =
                                if t2.IsNone then
                                    let ft3 = (fun t -> t*vCP) >> isSvalid
                                    let t3 = SearchDescendent(1000,0.0,1.0)(ft3)
                                    t3.Value*vCP
                                else 
                                    (1.0-t2.Value)*vCP+t2.Value*Vn
                                
                            
                            if mf(vDL)< mf(t1*Vn) then 
                                vDL
                            else
                                t1*Vn
                        
                        
                |None ->
                    failwith ""
            
            let vx = v.SubMatrix(0,x.Count,0,1)
            let vs =  S * v.SubMatrix(x.Count,s.Count,0,1)
            
            Matrix<float>.Build
               .DenseOfColumnArrays([| 
                                    Array.concat([vx.AsColumnMajorArray();
                                                  vs.AsColumnMajorArray()]) 
                                              |])
                            

        let LagrangeCalculator(s:Vector<float>,Ah:Matrix<float>,Ag:Matrix<float>,gf:Vector<float>,mu:float)=
                    
            let S = Matrix<float>.Build.DenseOfDiagonalVector(s)
            let Ak = Matrix<float>.Build.Dense(Ah.RowCount+S.RowCount,Ah.ColumnCount+Ag.ColumnCount)
            Ak.SetSubMatrix(0,0,Ah)
            Ak.SetSubMatrix(0,Ah.ColumnCount,Ag)
            Ak.SetSubMatrix(Ah.RowCount,Ah.ColumnCount,S)
            
            let ue = Vector<float>.Build.Dense(s.Count)+mu
            
            let vg = Vector<float>.Build.Dense(gf.Count+s.Count)
            vg.SetSlice(Some(0),None,-gf)
            vg.SetSlice(Some(gf.Count),None,ue)
                        
            
            let l = (Ak.Transpose()*Ak).Inverse()*Ak*vg
            
            let lh = l.SubVector(0,Ah.RowCount)
            let lg = l.SubVector(Ah.RowCount,Ag.RowCount)
            
            //Todo: process lg to modify the Hessian_s diagonal matrix
                        
            l
            
    type Nitro(feval:ScalarFunction,feq:VectorFunction option,fineq:VectorFunction option)=
        
                                 
        inherit IOptimizer(feval)
            override self.feq = feq
            override self.fineq = fineq
            override self.name="NITRO"
            
            override self.run(x0:Vector<float>)=
                let eta = 1e-8
                let tau = 0.995
                let theta = 0.2
                let zeta = 0.8
                let eu = 0.1
                let etol = 1e-7
                
                let mu = 0.1
                let delta0 = 1.0
                let nu0 =1.0
                
                let xk = x0.Clone()
                //let sk = Vector<float>.Build.DenseOfArray([|1.0,1.0,1.0|])
                
                
                let k = 0
                
//                while NITRO.E(xk,sk,0.0) > etol do  
//                    while NITRO.E(xk,sk,mu) > eu do
//                        //compute normal step
//                      ()
                                        
                
                failwith ""
                        
        