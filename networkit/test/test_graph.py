#!/usr/bin/env python3
import unittest
import random
import networkit as nk
import pickle

class TestGraph(unittest.TestCase):
	def getSmallGraph(self, weighted=False, directed=False):
		G = nk.Graph(4, weighted, directed)
		G.addEdge(0, 1, 1.0)
		G.addEdge(0, 2, 2.0)
		G.addEdge(3, 1, 4.0)
		G.addEdge(3, 2, 5.0)
		G.addEdge(1, 2, 3.0)

		return G

	def testNodeIterator(self):
		nk.setSeed(42, False)

		g = self.getSmallGraph()

		def doTest(g):
			nodes = []
			g.forNodes(lambda u: nodes.append(u))

			i = 0
			for u in g.iterNodes():
				self.assertEqual(u, nodes[i])
				i += 1

		doTest(g)
		g.removeNode(nk.graphtools.randomNode(g))
		g.removeNode(nk.graphtools.randomNode(g))
		doTest(g)

	def testEdgeIterator(self):
		for weighted in [True, False]:
			for directed in [True, False]:
				g = self.getSmallGraph(weighted, directed)
				for u, v in g.iterEdges():
					self.assertTrue(g.hasEdge(u, v))
				for u, v, w in g.iterEdgesWeights():
					self.assertTrue(g.hasEdge(u, v))
					self.assertEqual(g.weight(u, v), w)

	def testRemoveAllEdges(self):
		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)
				G.removeAllEdges()
				self.assertEqual(G.numberOfEdges(), 0)
				G.forNodePairs(lambda u, v: self.assertFalse(G.hasEdge(u, v)))

	def testRemoveSelfLoops(self):
		for directed in [True, False]:
			for weighted in [True, False]:
				g  = self.getSmallGraph(weighted, directed)
				for i in range(10):
					u = nk.graphtools.randomNode(g)
					g.addEdge(u, u)

				nSelfLoops = g.numberOfSelfLoops()
				nEdges = g.numberOfEdges()

				g.removeSelfLoops()

				self.assertEqual(nEdges - nSelfLoops, g.numberOfEdges())
				self.assertEqual(g.numberOfSelfLoops(), 0)

				g.forNodes(lambda u: self.assertFalse(g.hasEdge(u, u)))

	def testRemoveMultiEdges(self):
		nMultiedges = 5
		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)
				nEdges = G.numberOfEdges()

				for _ in range(nMultiedges):
					u, v = nk.graphtools.randomEdge(G)
					G.addEdge(u, v)

				G.removeMultiEdges()
				self.assertEqual(G.numberOfEdges(), nEdges)

				for _ in range(nEdges):
					u, v = nk.graphtools.randomEdge(G)
					G.removeEdge(u, v)
					self.assertFalse(G.hasEdge(u, v))
				self.assertEqual(G.numberOfEdges(), 0)


	def testNeighbors(self):
		# Directed
		G = nk.Graph(4, False, True)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)
		G.addEdge(1, 2)

		def getNeighbors(u):
			neighbors = []
			for v in G.iterNeighbors(u):
				neighbors.append(v)
			return sorted(neighbors)

		def getInNeighbors(u):
			inNeighbors = []
			for v in G.iterInNeighbors(u):
				inNeighbors.append(v)
			return sorted(inNeighbors)

		self.assertListEqual(getNeighbors(0), [1, 2])
		self.assertListEqual(getNeighbors(1), [2])
		self.assertListEqual(getNeighbors(2), [])
		self.assertListEqual(getNeighbors(3), [1, 2])

		self.assertListEqual(getInNeighbors(0), [])
		self.assertListEqual(getInNeighbors(1), [0, 3])
		self.assertListEqual(getInNeighbors(2), [0, 1, 3])
		self.assertListEqual(getInNeighbors(3), [])

		# Undirected
		G = nk.Graph(4, False, False)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)
		G.addEdge(1, 2)

		self.assertEqual(getNeighbors(0), [1, 2])
		self.assertEqual(getNeighbors(1), [0, 2, 3])
		self.assertEqual(getNeighbors(2), [0, 1, 3])
		self.assertEqual(getNeighbors(3), [1, 2])

	def testRandomEdgesReproducibility(self):
		numSamples = 10
		numSeeds = 3
		numRepeats = 10

		for directed in [False, True]:
			G = self.getSmallGraph(False, directed)

			results = [[] for i in range(numSeeds)]
			for repeats in range(numRepeats):
				for seed in range(numSeeds):
					nk.setSeed(seed, False)
					results[seed].append(nk.graphtools.randomEdges(G, numSamples))

			# assert results are different for different seeds
			for seed in range(1, numSeeds):
				self.assertNotEqual(results[0][0], results[seed][0])

			# assert results are consistent for same seeds
			for repeats in results:
				for x in repeats[1:]:
					self.assertEqual(repeats[0], x)

	def testIterator(self):
		# Undirected
		G = nk.Graph(4, False, False)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(1, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)

		def nbrFunc(u, v, weight, edgeId):
			forEdgesNbrs.append(v)

		# Iterate through neighbours of node using iterNeighbors
		def nodeIter(node):
			nodeList = []
			nbrs = G.iterNeighbors(node)
			try:
				while nbrs is not None:
					nodeList.append(next(nbrs))
			except StopIteration:
				pass
			return nodeList

		for node in range(G.upperNodeIdBound()):
			forEdgesNbrs = []
			G.forEdgesOf(node, nbrFunc)
			nodeNbrs = nodeIter(node)
			self.assertEqual(sorted(forEdgesNbrs), sorted(nodeNbrs))

		# Directed
		G = nk.Graph(4, False, True)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)
		G.addEdge(1, 2)

		# Iterate through neighbours of node using iterNeighbors
		def nodeInIter(node):
			nodeList = []
			nbrs = G.iterInNeighbors(node)
			try:
				while nbrs is not None:
					nodeList.append(next(nbrs))
			except StopIteration:
				pass
			return nodeList

		for node in range(G.upperNodeIdBound()):
			forEdgesNbrs = []
			G.forInEdgesOf(node, nbrFunc)
			nodeInNbrs = nodeInIter(node)
			self.assertEqual(sorted(forEdgesNbrs), sorted(nodeInNbrs))

	def testAddNodes(self):
		G = nk.Graph(0)
		G.addNodes(10)
		self.assertEqual(G.numberOfNodes(), 10)
		for u in range(G.numberOfNodes()):
			self.assertTrue(G.hasNode(u))

	def testAddEdgeWithMissingNodes(self):
		G = nk.Graph(0)

		with self.assertRaises(RuntimeError):
			G.addEdge(0, 1)

		self.assertEqual(G.numberOfNodes(), 0)
		self.assertEqual(G.numberOfEdges(), 0)

		G.addEdge(0, 2, addMissing = True)

		self.assertEqual(G.numberOfNodes(), 3)
		self.assertEqual(G.numberOfEdges(), 1)

		G.removeNode(1)

		self.assertEqual(G.numberOfNodes(), 2)
		self.assertEqual(G.numberOfEdges(), 1)

		G.addEdge(0, 1, addMissing = True)

		self.assertEqual(G.numberOfNodes(), 3)
		self.assertEqual(G.numberOfEdges(), 2)

		G.removeNode(2)

		G.addEdge(1, 2, addMissing = True)

		self.assertEqual(G.numberOfNodes(), 3)
		self.assertEqual(G.numberOfEdges(), 2)

	def testGraphPickling(self):
		G = nk.Graph(2)
		G.addEdge(0,1)
		pickle.dumps(G)

if __name__ == "__main__":
	unittest.main()
