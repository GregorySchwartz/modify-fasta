-- Diversity module.
-- By G.W. Schwartz
--
-- Collection of functions pertaining to finding the diversity of samples.

module Diversity where

-- Built-in
import Data.List

-- Takes two strings, returns Hamming distance
hamming :: String -> String -> Int
hamming xs ys = length $ filter not $ zipWith (==) xs ys

-- Returns the diversity of a list of things
diversity :: (Ord b) => Double -> [b] -> Double
diversity order sample
    | length sample == 0 = 0
    | order == 1         = exp . h $ speciesList
    | otherwise          = (sum . map ((** order) . p_i) $ speciesList) ** pow
  where
    pow          = 1 / (1 - order)
    h            = negate . sum . map (\x -> (p_i x) * (log (p_i x)))
    p_i x        = ((fromIntegral . length $ x) :: Double) /
                   ((fromIntegral . length $ sample) :: Double)
    speciesList  = group . sort $ sample

-- Calculates the binary coefficient
choose :: (Integral a) => a -> a -> a
choose _ 0 = 1
choose 0 _ = 0
choose n k = choose (n - 1) (k - 1) * n `div` k

-- Returns the rarefaction curve for each position in a list
rarefactionCurve :: (Eq a, Ord a) => [a] -> [Double]
rarefactionCurve xs = map rarefact [1..n_total]
  where
    rarefact n
        | n == 0       = 0
        | n == 1       = 1
        | n == n_total = k
        | otherwise    = k - ((1 / (fromIntegral (choose n_total n))) * inner n)
    inner n = fromIntegral                              .
              sum                                       .
              map (\g -> choose (n_total - length g) n) $
              grouped
    n_total = length xs
    k       = genericLength grouped
    grouped = group . sort $ xs

-- Calculates the percent of the curve that is above 95% of height of the curve
rarefactionViable :: [Double] -> Double
rarefactionViable xs = (genericLength valid / genericLength xs) * 100
  where
    valid = dropWhile (< (0.95 * last xs)) xs

