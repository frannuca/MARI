namespace MARIA.RBM
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.Random
type CRBM(Nv:int,Nh:int)=
       
    let mutable W = Matrix<float>.Build.Random(Nv+1,Nh+1,new MathNet.Numerics.Distributions.ContinuousUniform(0.0,1.0))/10.0
    let mutable dW = Matrix<float>.Build.Random(Nv+1,Nh+1,new MathNet.Numerics.Distributions.ContinuousUniform(0.0,1.0))/1000.0
  
    let mutable Av = Vector<float>.Build.DenseOfArray  [| for i in 0 .. Nv do yield 0.05|] //Exponential factor for sigmoids in visible layer
    let mutable Ah = Vector<float>.Build.DenseOfArray  [| for i in 0 .. Nh do yield 1.0|] //Exponential factor for sigmoids in hidden layer
    
    let dAh = Vector<float>.Build.DenseOfArray(Array.zeroCreate(Nh+1))

    let sigmoid = (new activation.Sigmoid(-5.0,5.0)) :> activation.IActivation<float>

    //Optimization parameters:
    let sigma = 0.2
    let epsW = 0.005
    let epsA = 0.005
    let cost = 0.0
    let momentum = 0.0

    let normd = new MathNet.Numerics.Distributions.Normal(0.0,sigma,new System.Random(42))
       
    let updateHidden(sv:Vector<float>)=      
        let e = Vector<float>.Build.Random(Nh+1,normd)
        let auxs = W.Transpose()*sv + e
        let sr = sigmoid.vec_eval(Ah.ToArray()) auxs 
        let n = sr.Count
        sr.[ n - 1] <- 1.0
        sr

    let updateVisible(sh:Vector<float>)=      
        let e = Vector<float>.Build.Random(Nv+1,normd)
        let auxs = W*sh + e
        let sr = sigmoid.vec_eval(Av.ToArray()) auxs 
        let n = sr.Count
        sr.[ n - 1] <- 1.0
        sr
        
    
    member self.learn(data:Matrix<float>)(niter:int,k:int)=
        let mutable err = Array.zeroCreate<float>(niter)
        let mutable E = Array.zeroCreate<float>(niter)

        let nsamples = float(data.ColumnCount)
        for i in 0 .. niter-1 do
            
            let mutable wpos = Matrix<float>.Build.Dense(Nv+1,Nh+1)
            let mutable wneg = Matrix<float>.Build.Dense(Nv+1,Nh+1)
            let mutable apos = Vector<float>.Build.Dense(Nh+1)
            let mutable aneg = Vector<float>.Build.Dense(Nh+1)
            let mutable sh= Vector<float>.Build.Dense(Nh+1)
            for n in 0 .. int(nsamples)-1 do
                let xin = Vector<float>.Build.DenseOfArray(Array.concat [data.Column(n).ToArray();[|1.0|]])
                
                sh <- updateHidden(xin)
                wpos <- wpos + xin.OuterProduct(sh)
                apos <- apos + sh.PointwiseMultiply(sh)

                let mutable sv = updateVisible(sh)
                sh <- updateHidden(sv)

                for i in 0 .. k do
                    sv <- updateVisible(sh)
                    sh <- updateHidden(sv)
                
                let tmp = sv.OuterProduct(sh)
                wneg <- wneg + tmp
                aneg <- aneg + sh.PointwiseMultiply(sh)

                let delta = sv - xin
                err.[i] <- err.[i]+ delta.L2Norm()
            
            printfn " iter=%i, ERROR=%f" i err.[i]
            dW <- dW * momentum + epsW* (wpos-wneg)/nsamples - cost*W
            W <- W + dW
            Ah <- Ah + epsA* (apos-aneg)/(nsamples* sh.PointwiseMultiply(sh))
    
    member self.reconstruct(nPoints:int,nSteps:int)=
        let Z = Matrix<float>.Build.Dense(Nv,nPoints)

        for i in 0 .. nPoints-1 do
            let xr = Vector<float>.Build.Random(Nv)
            let x = Vector<float>.Build.DenseOfArray(Array.concat [xr.ToArray();[|1.0|]])
            
            let mutable sh = updateHidden(x)
            let mutable sv = updateVisible(sh)
            for i in 0 .. nPoints-2 do
                sh <- updateHidden(sv)
                sv <- updateVisible(sh)
            
            Z.SetColumn(i,sv.SubVector(0,sv.Count-1))    
        Z
    
