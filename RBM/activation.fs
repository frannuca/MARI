namespace MARIA.RBM
open MathNet.Numerics.LinearAlgebra

module activation=
    type IActivation<'T>=
        abstract eval:'T -> float -> float
        abstract vec_eval: 'T array -> Vector<float> -> Vector<float>
    

    type Sigmoid(lo:float,hi:float)= 
        let lo_ = lo
        let hi_ = hi
        let feval a x = lo_ + (hi_ - lo_) / (1.0 + exp(-a * x))
        interface IActivation<float> with
            member self.eval (a:float)(x:float)=
                feval a x

            member self.vec_eval(A:array<float>)(x:Vector<float>)=               
               Vector.Build.DenseOfArray( A |> Array.mapi(fun i a -> feval(a)(x.[i]))) 
                
            

