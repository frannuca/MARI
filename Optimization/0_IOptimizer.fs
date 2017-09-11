namespace MARI.Optimization

open MathNet.Numerics.LinearAlgebra




    type ScalarFunction = Vector<float> -> float
    type VectorFunction = Vector<float> -> Vector<float>
    type MatrixFcuntion = Vector<float> -> Matrix<float>
    
    type FContraint = VectorFunction option
         
    
    type Result={nevals:int;minfeval:float;xsolution:Vector<float>; ftol:float;xtol:float}
    
    [<AbstractClass>]
    type IOptimizer(f:ScalarFunction) =
        abstract member feq:FContraint
        abstract member fineq:FContraint        
        
        abstract member name:string
        
        member  self.feval = f
        
        default self.feq = None
        default self.fineq = None
        
        
        abstract member run: Vector<float> -> Result
        
        
                
        
    
    
        