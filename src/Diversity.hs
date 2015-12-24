-- Diversity module.
-- By G.W. Schwartz
--
-- Collection of functions pertaining to finding the diversity of samples.

module Diversity where

-- Built-in
import Data.List
import qualified Data.Text as T

-- | Takes two strings, returns Hamming distance
hamming :: T.Text -> T.Text -> Int
hamming xs ys = length $ filter (not . uncurry (==)) $ T.zip xs ys
