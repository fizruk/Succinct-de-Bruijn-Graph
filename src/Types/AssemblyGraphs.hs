-- |
-- For asymptotic analysis below we use the following variables:
--
-- * p - base of the graph
-- * n - sdf
-- * m - asdf
module Types.AssemblyGraphs where

import           Control.Arrow ((&&&))
import           Data.Array
import           Data.BitArray
import           Data.Function (on)
import qualified Data.Set      as Set
import           Prelude       hiding (seq)
import           Types.DNA

-- $setup
-- >>> :set -XOverloadedStrings

type Edge = Int

type Node = Int

-- | Please write full sentences.
type Base = Int

-- | Please write full sentences.
data DeBruijnGraph = DeBruijnGraph
  { p           :: Base             -- ^ Please write full sentences.
  , bitArr      :: BitArray         -- ^ Please write full sentences.
  , includes    :: [(Int, Bool)]
  , occurrences :: Array Edge Int
  }

instance Show DeBruijnGraph where
  show deBruijnGraph = show multiplicityList
    where
      multiplicityList =
        map
          (\(ind, value) -> (numberToSequence (p deBruijnGraph) ind, value))
          ((assocs . occurrences) deBruijnGraph)

-- | Graphs are equal when
--
-- * base ...
-- * bitarray ...
-- * occurrences ...
instance Eq DeBruijnGraph where
--  g1 == g2 = and
--    [ p g1           == p g2
--    , bitArr g1      == bitArr g2
--    , occurrences g1 == occurrences g2
--    ]
  (==) = (==) `on` (p &&& bitArr &&& occurrences)

-- * Construction

-- | Empty de Bruijn graph.
--
-- >>> emptyDeBruijn 1
-- [(A,0),(C,0),(G,0),(T,0)]
emptyDeBruijn :: Base -> DeBruijnGraph
emptyDeBruijn base =
  DeBruijnGraph
    { p = base
    , bitArr = bitArray (0, 4 ^ base - 1) []
    , includes = []
    , occurrences = array (0, 4 ^ base - 1) [(i, 0) | i <- [0 .. 4 ^ base - 1]]
    }

-- | /O(n)/.
-- Add DNA sequence to a de Bruijn graph.
insertSequence :: DNASequence -> DeBruijnGraph -> DeBruijnGraph
insertSequence seq deBruijnGraph
  | length (show seq) > p deBruijnGraph =
    insertSequences (splitByN (p deBruijnGraph) seq) deBruijnGraph
  | otherwise =
    deBruijnGraph
      {bitArr = newBitArray, includes = newInserted, occurrences = newArr}
  where
    setNumber = sequenceToNumber seq
    newInserted = (setNumber, True) : includes deBruijnGraph
    newBitArray = bitArray (bitArrayBounds (bitArr deBruijnGraph)) newInserted
    occures = (occurrences deBruijnGraph ! setNumber) + 1
    newArr = occurrences deBruijnGraph // [(setNumber, occures)]

insertSequences :: [DNASequence] -> DeBruijnGraph -> DeBruijnGraph
insertSequences [] deBruijnGraph = deBruijnGraph
insertSequences (seq:seqs) deBruijnGraph = insertSequences seqs newDeBruijnGraph
  where
    newDeBruijnGraph = insertSequence seq deBruijnGraph

fromDNASequences :: Base -> [DNASequence] -> DeBruijnGraph
fromDNASequences base seqs = insertSequences seqs (emptyDeBruijn base)

diffSequence :: DeBruijnGraph -> DNASequence -> DeBruijnGraph
diffSequence deBruijnGraph seq' =
  deBruijnGraph
    {bitArr = newBitArr, includes = newIncludes, occurrences = newOccurrences}
  where
    edge = sequenceToNumber seq'
    occurrence = occurrences deBruijnGraph ! edge
    newOccurrences = occurrences deBruijnGraph // [(edge, occurrence - 1)]
    newBitArr = bitArray (0, 4 ^ p deBruijnGraph - 1) newIncludes
    newIncludes =
      if newOccurrences ! edge == 0
        then filter (\(ind, _) -> ind /= edge) (includes deBruijnGraph)
        else includes deBruijnGraph

(///) :: DeBruijnGraph -> DNASequence -> DeBruijnGraph
(///) = diffSequence

eulerBackTracking ::
     DeBruijnGraph -> Edge -> [Edge] -> [DNASequence] -> [DNASequence]
eulerBackTracking deBruijnGraph (-1) (newCurrent:xs) path =
  eulerBackTracking deBruijnGraph newCurrent xs path
eulerBackTracking deBruijnGraph current [] path
  | null successors = newPath
  | otherwise = eulerPath deBruijnGraph successor [] newPath
  where
    successors =
      (successorEdges (bitArr deBruijnGraph) .
       getFromNode . numberToSequence (p deBruijnGraph))
        current
    successor = head successors
    newPath = numberToSequence (p deBruijnGraph) current : path
eulerBackTracking deBruijnGraph current visited@(newCurrent:xs) path
  | null successors = eulerBackTracking deBruijnGraph newCurrent xs newPath
  | otherwise = eulerPath deBruijnGraph successor visited newPath
  where
    successors =
      (successorEdges (bitArr deBruijnGraph) .
       getFromNode . numberToSequence (p deBruijnGraph))
        current
    successor = head successors
    newPath = numberToSequence (p deBruijnGraph) current : path

-- |
eulerPath
  :: DeBruijnGraph    -- ^
  -> Edge             -- ^
  -> [Edge]
  -> [DNASequence]
  -> [DNASequence]  -- TODO: change to Maybe [DNASequence]
eulerPath deBruijnGraph current visited path =
  case successors of
    (newCurrent:_) -> eulerPath newDeBruijnGraph newCurrent newVisited path
    _              -> eulerBackTracking newDeBruijnGraph (-1) newVisited path
  where
    successors =
      (successorEdges (bitArr newDeBruijnGraph) .
       getToNode . numberToSequence (p newDeBruijnGraph))
        current
    newDeBruijnGraph =
      deBruijnGraph /// numberToSequence (p deBruijnGraph) current
    newVisited = current : visited

selectStartNode :: [(Node, (Int, Int))] -> Node
selectStartNode l = startNode
  where
    calculatedL =
      map (\(ind, (inEdges, outEdges)) -> (ind, outEdges - inEdges)) l
    filteredL = filter (\(_, diff) -> diff == 1) calculatedL
    startNode =
      if null filteredL
        then fst $ head (filter (\(_, (_, out)) -> out > 0) l)
        else fst $ head filteredL

selectEndNode :: [(Node, (Int, Int))] -> Node
selectEndNode l = endNode
  where
    calculatedL =
      map (\(ind, (inEdges, outEdges)) -> (ind, inEdges - outEdges)) l
    filteredL = filter (\(_, diff) -> diff == 1) calculatedL
    endNode =
      if null filteredL
        then fst $ head (filter (\(_, (in', _)) -> in' > 0) l)
        else fst $ head filteredL

calculateInEdges :: Array Edge Int -> Node -> Base -> Int
calculateInEdges arr node base =
  arr ! node + arr ! (node + 4 ^ base) + arr ! (node + 4 ^ base * 2) +
  arr ! (node + 4 ^ base * 3)

calculateOutEdges :: Array Edge Int -> Node -> Int
calculateOutEdges arr node =
  arr ! node + arr ! (node + 1) + arr ! (node + 2) + arr ! (node + 3)

selectNodes :: Array Edge Int -> Base -> (Node, Node)
selectNodes arr base =
  (selectStartNode $ assocs inOutArray, selectEndNode $ assocs inOutArray)
  where
    inOutArray = array (0, arrSize) [inOut i | i <- [0 .. arrSize]]
    arrSize = (4 ^ (base - 1)) - 1
    inOut i =
      (i, (calculateInEdges arr i (base - 1), calculateOutEdges arr (4 * i)))

selectStartEdge :: DeBruijnGraph -> Edge
selectStartEdge deBruijnGraph =
  head $
  successorEdges (bitArr deBruijnGraph) $
  fst (selectNodes (occurrences deBruijnGraph) (p deBruijnGraph))

successorEdges :: BitArray -> Node -> [Edge]
successorEdges bitArr' node = s
  where
    ranks =
      filter
        (> rank bitArr' (4 * node - 1))
        [rank bitArr' (4 * node + n) | n <- [0 .. 3]]
    successors = map (select bitArr') ranks
    s = filter (>= 0) $ Set.toList $ Set.fromList successors

select :: BitArray -> Int -> Int
select bitArr' i = select' bitList i 0 - 1
  where
    select' [] _ ind          = ind
    select' (True:_) 0 ind    = ind + 1
    select' (False:_) 0 ind   = ind
    select' (True:xs) i' ind  = select' xs (i' - 1) (ind + 1)
    select' (False:xs) i' ind = select' xs i' (ind + 1)
    bitList = bits bitArr'

rank :: BitArray -> Int -> Int
rank bitArr' i = sum $ take (i + 1) bitList
  where
    bitList = bits01 bitArr'

assemblyDeBruijn :: DeBruijnGraph -> DNASequence
assemblyDeBruijn deBruijnGraph = foldl (+++) (DNASequence "") eulerPath'
  where
    eulerPath' = eulerPath deBruijnGraph startEdge [] []
    startEdge = selectStartEdge deBruijnGraph
