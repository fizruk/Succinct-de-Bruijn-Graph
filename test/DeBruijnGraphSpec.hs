module DeBruijnGraphSpec where

import           Test.Hspec
import           Types.AssemblyGraphs
import           Types.DNA

dnaSequences :: [DNASequence]
dnaSequences = map DNASequence ["TT", "TC", "CG", "GG", "GA", "AA", "AG"]

deBruijnGraph :: DeBruijnGraph
deBruijnGraph = fromDNASequences 2 dnaSequences

testBuildSimpleGraph :: IO ()
testBuildSimpleGraph = assembledSequence `shouldBe` DNASequence "ATGGCGTGCA"
  where
    assembledSequence = assemblyDeBruijn deBruijnGraph'
    dnaSequences' =
      map DNASequence ["CGTGCAA", "ATGGCGT", "CAATGGC", "GGCGTGC", "TGCAATG"]
    deBruijnGraph' = fromDNASequences 3 dnaSequences'

testSuccessor :: IO ()
testSuccessor = successors `shouldBe` map DNASequence ["TC", "TT"]
  where
    successors =
      map (numberToSequence 2) $
      successorEdges (bitArr deBruijnGraph) (sequenceToNumber (DNASequence "T"))

testSequenceToNumber :: IO ()
testSequenceToNumber = num `shouldBe` 10
  where
    num = sequenceToNumber (DNASequence "GG")

testSelect :: IO ()
testSelect = s `shouldBe` 0
  where
    s = select (bitArr deBruijnGraph) 1

testRank :: IO ()
testRank = s `shouldBe` 7
  where
    s = rank (bitArr deBruijnGraph) 40

testSelectStartEdge :: IO ()
testSelectStartEdge = numberToSequence 2 e `shouldBe` DNASequence "TC"
  where
    e = selectStartEdge deBruijnGraph

testSelectNodes :: IO ()
testSelectNodes = nodes `shouldBe` (3, 2)
  where
    nodes = selectNodes (occurrences deBruijnGraph) (p deBruijnGraph)

testToNode :: IO ()
testToNode = numberToSequence 1 n `shouldBe` DNASequence "C"
  where
    n = getToNode $ DNASequence "TC"

testEulerBackTracking :: IO ()
testEulerBackTracking = s `shouldBe` map DNASequence ["TT", "TC"]
  where
    s =
      eulerBackTracking
        graph
        ((sequenceToNumber . DNASequence) "TT")
        []
        [DNASequence "TC"]
    graph = emptyDeBruijn 2

testEulerPath :: IO ()
testEulerPath = s `shouldBe` map DNASequence ["TT"]
  where
    s = eulerPath graph startEdge [] []
    startEdge = (sequenceToNumber . DNASequence) "TT"
    graph = fromDNASequences 2 [DNASequence "TT"]

testCheckEdge :: IO ()
testCheckEdge = g `shouldBe` emptyDeBruijn 2
  where
    g = graph /// DNASequence "TT"
    graph = fromDNASequences 2 [DNASequence "TT"]

spec :: Spec
spec =
  describe "Tests for de Bruijn Graph" $ do
    it "Sequence to number" testSequenceToNumber
    it "Select" testSelect
    it "Rank" testRank
    it "Successor" testSuccessor
    it "Select Pivoting Nodes" testSelectNodes
    it "Select Start Edge" testSelectStartEdge
    it "Check edge" testCheckEdge
    it "Euler Back Tracking" testEulerBackTracking
    it "Euler Path" testEulerPath
    it "Simple cyclic genome" testBuildSimpleGraph
