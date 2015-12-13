import Data.Word
import Data.List(elemIndex)
import Data.Maybe
import Data.Complex
import qualified Data.Array.Repa as R
import Data.Array.Repa.IO.BMP
import System.Environment

type Number = Double
type ComplexNumber = Complex Number
type IterationCount = Integer
type PointData = (ComplexNumber, IterationCount)
type RootData = (Int, IterationCount)
type Bounds = (Number, Number, Number, Number)
type Step = Number
type Color = (Word8, Word8, Word8)

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
epsilon = 1e-12

newton :: Function -> Function -> ComplexNumber -> PointData
newton f f' z0 = newton' z1 z0 0
    where approximate z = z - (evaluate f z) / (evaluate f' z)
          z1 = approximate z0
          newton' current previous n
              | n > maxIterations = (current, n)
              | magnitude (current - previous) < epsilon = (current, n)
              | otherwise = newton' (approximate current) current (n + 1)

windowApply :: (ComplexNumber -> PointData) -> Bounds -> Step -> [PointData]
windowApply f (top, left, bottom, right) step = [f (x :+ y ) | y <- ys, x <- xs]
    where xs = [left, left + step .. right]
          ys = [top, top + step .. bottom]

approximateIndex :: ComplexNumber -> [ComplexNumber] -> Maybe Int
approximateIndex z zs = elemIndex True [magnitude (z - zi) < epsilon * 2 | zi <- zs]

closestIndex :: ComplexNumber -> [ComplexNumber] -> Int
closestIndex z zs = snd . minimum $ zip [magnitude (z - zi) | zi <- zs] [0..]

closeEnough :: ComplexNumber -> [ComplexNumber] -> Bool
closeEnough z zs = case approximateIndex z zs of
    Just _ -> True
    Nothing ->  False

roots :: [PointData] -> [ComplexNumber]
roots ps = foldr checkRoot [] zs
    where checkRoot z roots
              | z `closeEnough` roots = roots
              | otherwise = z:roots
          zs = fst . unzip $ ps

assignRoot :: [ComplexNumber] -> PointData -> RootData
assignRoot rs p@(z, i) = (closestIndex z rs, i)

colorRoot :: RootData -> Color
colorRoot (n, i) = (ra, ga, ba)
    where colors = [ (0, 0, 0)
                   , (255, 0, 0)
                   , (0, 255, 0)
                   , (0, 0, 255)
                   , (255, 255, 0)
                   , (255, 0, 255)
                   , (0, 255, 255)
                   , (255, 255, 255)
                   ]
          (r, g, b) = colors !! n
          ii = max i 1
          ra = fromIntegral . min 255  $ 4 * r `div` ii
          ga = fromIntegral . min 255  $ 4 * g `div` ii
          ba = fromIntegral . min 255  $ 4 * b `div` ii


colorRoots :: Function -> Bounds -> Step -> [Color]
colorRoots f bounds step = map colorRoot rds
    where ps = windowApply (newton f (derivative f)) bounds step
          rs = roots ps
          rds = map (assignRoot rs) ps

arrayify :: Int -> Int -> [Color] -> R.Array R.U R.DIM2 Color
arrayify width height =
    R.fromListUnboxed (R.Z R.:. height R.:. width :: R.DIM2)

f = Sum (Power 3 Z) (realConstant (-1))

main = do
    args <- getArgs
    let [fileName, stringTop, stringLeft, stringBottom, stringRight, stringStep] = args
        top = read stringTop
        left = read stringLeft
        bottom = read stringBottom
        right = read stringRight
        step = read stringStep
        width = floor ((right - left) / step) + 1
        height = floor ((bottom - top) / step) + 1
        colorList = colorRoots f (top, left, bottom, right) step
        colorArray = arrayify width height colorList
    writeImageToBMP fileName colorArray
