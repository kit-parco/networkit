""" This module handles community detection, i.e. the discovery of densely connected groups in networks."""

__author__ = "Christian Staudt"


from _NetworKit import Partition, Coverage, Modularity, CommunityDetector, PLP, LPDegreeOrdered, PLM, PartitionReader, PartitionWriter,\
	NodeStructuralRandMeasure, GraphStructuralRandMeasure, JaccardMeasure, NMIDistance, AdjustedRandMeasure,\
	StablePartitionNodes, IntrapartitionDensity, PartitionHubDominance, CoverHubDominance, PartitionFragmentation, IsolatedInterpartitionExpansion, IsolatedInterpartitionConductance,\
	EdgeListPartitionReader, GraphClusteringTools, ClusteringGenerator, PartitionIntersection, HubDominance, CoreDecomposition, CutClustering, ParallelPartitionCoarsening, \
	BinaryPartitionReader, BinaryPartitionWriter, BinaryEdgeListPartitionReader, BinaryEdgeListPartitionWriter, \
	CoverF1Similarity ,EgoSplitting, LPPotts, LouvainMapEquation, \
	PLMFactory, PLPFactory, LouvainMapEquationFactory, LPPottsFactory

# R.I.P.: The CNM (Clauset, Newman, Moore) community detection algorithm - it was always a bit slow, but it broke down in the end. Resurrect it from history (<= 3.4.1) if needed for experimental purposes.

# local imports
#from .properties import CoreDecomposition, overview
from . import graph
from . import stopwatch
from . import graphio
from .support import MissingDependencyError

# external imports
import os
import math
import random
import warnings
try:
	import tabulate
except ImportError:
	have_tabulate = False
else:
	have_tabulate = True
import tempfile
import subprocess

def detectCommunities(G, algo=None, inspect=True):
	""" Perform high-performance community detection on the graph.
		:param    G    the graph
		:param     algorithm    community detection algorithm instance
		:return communities (as type Partition)
		"""
	if algo is None:
		algo = PLM(G, refine=False)
	t = stopwatch.Timer()
	algo.run()
	zeta = algo.getPartition()
	t.stop()
	print("{0} detected communities in {1} [s]".format(algo.toString(), t.elapsed))
	if inspect:
		print ("solution properties:")
		inspectCommunities(zeta, G)
	return zeta

def inspectCommunities(zeta, G):
	""" Display information about communities
		:param    zeta    communities
		:param    G        graph
	"""
	if not have_tabulate:
		raise MissingDependencyError("tabulate")
	communitySizes = zeta.subsetSizes()
	mod = Modularity().getQuality(zeta, G)
	commProps = [
		["# communities", zeta.numberOfSubsets()],
		["min community size", min(communitySizes)],
		["max community size", max(communitySizes)],
		["avg. community size", sum(communitySizes) / len(communitySizes)],
		#["imbalance", zeta.getImbalance()],
		["modularity", mod],
	]
	print(tabulate.tabulate(commProps))


def communityGraph(G, zeta):
	""" Create a community graph, i.e. a graph in which one node represents a community and an edge represents the edges between communities, from a given graph and a community detection solution"""
	cg = ParallelPartitionCoarsening(G, zeta)
	cg.run()
	return cg.getCoarseGraph()


def evalCommunityDetection(algo, G):
	""" Evaluate a community detection algorithm """

	if not have_tabulate:
		raise MissingDependencyError("tabulate")
	t = stopwatch.Timer()
	algo.run()
	zeta = algo.getPartition()
	t.stop()
	results = [
		["time [s]", t.elapsed],
		["# communities", zeta.numberOfSubsets()],
		["modularity", Modularity().getQuality(zeta, G)]
	]
	print(tabulate.tabulate(results))


def readCommunities(path, format="default"):
	""" Read a partition into communities from a file"""
	readers =  {"default": PartitionReader(),
		"edgelist-t1": EdgeListPartitionReader(1, '\t'),
		"edgelist-t0": EdgeListPartitionReader(0, '\t'),
		"edgelist-s1": EdgeListPartitionReader(1, ' '),
		"edgelist-s0": EdgeListPartitionReader(0, ' '),
		}
	# get reader
	try:
		reader = readers[format]#(**kwargs)
	except KeyError:
		raise Exception("unrecognized format: {0}".format(format))

	# get proper file path
	if ("~" in path):
		path = os.path.expanduser(path)
		print("path expanded to: {0}".format(path))
	# check if file path leads to a valid file
	if not os.path.isfile(path):
		raise IOError("{0} is not a file".format(path))
	else:
		with open(path, "r") as file:    # catch a wrong path before it crashes the interpreteri
			print("read communities from: {0}".format(path))
			communities = reader.read(path)
			return communities

	return None


def writeCommunities(communities, path):
	""" Write a partition into communities to a file"""
	PartitionWriter().write(communities, path)
	print("wrote communities to: {0}".format(path))


def compareCommunities(G, zeta1, zeta2):
	""" Compare the partitions with respect to several (dis)similarity measures"""
	raise NotImplementedError("TODO:")


def kCoreCommunityDetection(G, k, algo=None, inspect=True):
	""" Perform community detection on the k-core of the graph, which possibly
		reduces computation time and enhances the result.
		:param    G    the graph (may not contain self-loops)
		:param		k 	k as in k-core
		:param     algorithm    community detection algorithm instance
		:return communities (as type Partition)
		"""
	coreDec = CoreDecomposition(G)
	coreDec.run()

	cores = coreDec.cores()
	try:
		kCore = cores[k]
	except IndexError:
		raise Error("There is no core for the specified k")

	C = graph.Subgraph().fromNodes(G, kCore)	# FIXME: node indices are not preserved

	#properties.overview(C)

	return detectCommunities(C, algo, inspect)


def mesoscopicResponseFunction(G, samples=100):
	"""
	"""
	raise NotImplementedError("work in progress")
	m = G.numberOfEdges()
	gammaRangeLow = [math.e**x for x in range(-10, 0)]
	gammaRangeHigh = [math.e**x for x in range(0, math.ceil(math.log(2*m)))]
	gammaRange = gammaRangeLow + gammaRangeHigh
	print(gammaRange)
	nCom = []

	for gamma in gammaRange:
		communityDetector = PLM(G, gamma=gamma)
		communityDetector.run()
		communities = communityDetector.getPartition()
		nCom.append(communities.numberOfSubsets())

	return (gammaRange, nCom)


class InfomapAdapter:

	infomapPath = None

	def __init__(self, G):
		self.G = G

	@classmethod
	def setPath(cls, infomapPath):
		cls.infomapPath = infomapPath

	def run(self):
		if not self.infomapPath:
			raise Exception("set path to infomap binary with 'setPath' class method")
		with tempfile.TemporaryDirectory() as tempdir:
			print("temporary file directory: ", tempdir)
			graph_filename = os.path.join(tempdir, "network.txt")
			graphio.writeGraph(self.G, graph_filename, fileformat=graphio.Format.EdgeListSpaceZero)
			subprocess.call([self.infomapPath, "-s", str(random.randint(-2**31, 2**31)), "-2", "-z", "--clu", graph_filename, tempdir])
			self.result = readCommunities(os.path.join(tempdir, "network.clu"), format="edgelist-s0")
			while self.result.numberOfElements() < self.G.upperNodeIdBound():
				self.result.toSingleton(result.extend())
		return self

	def getPartition(self):
		return self.result
