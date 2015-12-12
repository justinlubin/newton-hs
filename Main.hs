import Data.List(elemIndex)
import Data.Complex

type Number = Double
type ComplexNumber = Complex Number
type IterationCount = Integer
type DataPoint = (ComplexNumber, IterationCount)
type Bounds = (Number, Number, Number, Number)

e :: ComplexNumber
e = (exp 1) :+ 0

data Function = Constant ComplexNumber
              | Z
              | Sum Function Function
              | Product Function Function
              | Quotient Function Function
              | Power ComplexNumber Function
              | Sin Function
              | Cos Function
              | Exponential ComplexNumber Function
              | Log ComplexNumber Function
              | Undefined
              | Indeterminate
              deriving (Show)

realConstant :: Number -> Function
realConstant n = Constant (n :+ 0)

derivative :: Function -> Function
derivative (Constant _) = realConstant 0
derivative Z = realConstant 1
derivative (Sum u v) = Sum (derivative u) (derivative v)
derivative (Product u v) = Sum (Product (derivative u) v) (Product u (derivative v))
derivative (Quotient u v) = (Quotient (Sum (Product v (derivative u)) (Product (Product (realConstant (-1)) u) (derivative v))) (Power 2 v))
derivative (Power p f) = Product (Product (Constant p) (Power (p - 1) f)) (derivative f)
derivative (Sin f) = Product (Cos f) (derivative f)
derivative (Cos f) = Product (Product (realConstant (-1)) (Sin f)) (derivative f)
derivative (Exponential b f) = Product (Product (Exponential b f) (Log e (Constant b))) (derivative f)
derivative (Log b f) = Quotient (derivative f) (Product f (Log e (Constant b)))
derivative Undefined = Undefined
derivative Indeterminate = Indeterminate

evaluate :: Function -> ComplexNumber -> ComplexNumber
evaluate (Constant a) z = a
evaluate Z z = z
evaluate (Sum u v) z = (evaluate u z) + (evaluate v z)
evaluate (Product u v) z = (evaluate u z) * (evaluate v z)
evaluate (Quotient u v) z = (evaluate u z) / (evaluate v z)
evaluate (Power p f) z = (evaluate f z) ** p
evaluate (Sin f) z = sin (evaluate f z)
evaluate (Cos f) z = cos (evaluate f z)
evaluate (Exponential b f) z = b ** (evaluate f z)
evaluate (Log b f) z = logBase b (evaluate f z)
evaluate Undefined z = undefined
evaluate Indeterminate z = undefined

maxIterations :: IterationCount
maxIterations = 100

epsilon :: Number
epsilon = 1e-9

newton :: Function -> ComplexNumber -> DataPoint
newton f z0 = newton' f z1 z0 0
    where f' = derivative f
          approximate f z = z - (evaluate f z) / (evaluate f' z)
          z1 = approximate f z0
          newton' f current previous n
              | n > maxIterations = (current, n)
              | magnitude (current - previous) < epsilon = (current, n)
              | otherwise = newton' f (approximate f current) current (n + 1)

windowApply :: (ComplexNumber -> DataPoint) -> Bounds -> Number -> [[DataPoint]]
windowApply f (top, left, bottom, right) step = [[f (x :+ y ) | x <- xs] | y <- ys]
    where xs = [left, left + step .. right]
          ys = [top, top + step .. bottom]

approximateIndex :: ComplexNumber -> [ComplexNumber] -> Maybe Int
approximateIndex z zs = elemIndex True [magnitude (z - zi) < 1 | zi <- zs]

closeEnough :: ComplexNumber -> [ComplexNumber] -> Bool
closeEnough z zs = case approximateIndex z zs of
    Just _ -> True
    Nothing -> False

roots :: [[DataPoint]] -> [ComplexNumber]
roots field = foldr checkRoot [] flat
    where checkRoot z roots
              | z `closeEnough` roots = z:roots
              | otherwise = roots
          flat = fst . unzip . concat $ field

--assignRoots :: Function -> Bounds -> Number -> [[(Int, IterationCount)]]
assignRoots f bounds step = assignRoots' [] windowApplication
    where windowApplication = windowApply (newton f) bounds step
          assignRoot roots z = case ai of
              Just i -> (i, roots)
              Nothing -> (length roots, roots ++ [z])
              where ai = approximateIndex z roots
          assignRoots' roots [[]] = []
          assignRoots' roots ([]:rows) = assignRoots' roots rows
          assignRoots' roots ((tuple:tuples):rows) = (fst ar, snd tuple,"      ") : assignRoots' (snd ar) (tuples:rows)
              where ar = assignRoot roots (fst tuple)

hi = Sum (Power 3 Z) (realConstant (-1))

main = print $ roots $ windowApply (newton hi) (-10, -10, 10, 10) 0.1
