# distutils: language=c++

from libc.stdint cimport uint64_t

from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.map cimport map

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

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node
ctypedef double edgeweight

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition, _Cover, Cover
from .graphio import PartitionReader, PartitionWriter, EdgeListPartitionReader, BinaryPartitionReader, BinaryPartitionWriter, BinaryEdgeListPartitionReader, BinaryEdgeListPartitionWriter

from . import graph
from .centrality import CoreDecomposition
from .coarsening import ParallelPartitionCoarsening
from . import stopwatch
from . import graphio
from .support import MissingDependencyError

cdef extern from "<algorithm>" namespace "std":
	pair[_Graph, vector[node]] move(pair[_Graph, vector[node]]) nogil

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

none = _none

cdef extern from "<networkit/community/CommunityDetectionAlgorithm.hpp>":

	cdef cppclass _CommunityDetectionAlgorithm "NetworKit::CommunityDetectionAlgorithm"(_Algorithm):
		_CommunityDetectionAlgorithm(const _Graph &_G)
		_Partition getPartition() except +

cdef class CommunityDetector(Algorithm):
	""" Abstract base class for static community detection algorithms """

	cdef Graph _G
	def __init__(self, *args, **namedargs):
		if type(self) == CommunityDetector:
			raise RuntimeError("Error, you may not use CommunityDetector directly, use a sub-class instead")

	def getPartition(self):
		"""  Returns a partition of the clustering.

		Returns:
		--------
		networkit.Partition:
			A Partition of the clustering.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return Partition().setThis((<_CommunityDetectionAlgorithm*>(self._this)).getPartition())

# Fused type for methods that accept both a partition and a cover
ctypedef fused PartitionCover:
	Partition
	Cover

cdef extern from "<networkit/community/ClusteringGenerator.hpp>":

	cdef cppclass _ClusteringGenerator "NetworKit::ClusteringGenerator":
		_ClusteringGenerator() except +
		_Partition makeSingletonClustering(_Graph G) except +
		_Partition makeOneClustering(_Graph G) except +
		_Partition makeRandomClustering(_Graph G, count k) except +
		_Partition makeContinuousBalancedClustering(_Graph G, count k) except +
		_Partition makeNoncontinuousBalancedClustering(_Graph G, count k) except +

cdef class ClusteringGenerator:
	""" Generators for various clusterings """
	cdef _ClusteringGenerator _this
	def makeSingletonClustering(self, Graph G):
		"""  Generate a clustering where each node has its own cluster

		Parameters:
		-----------
		G : networkit.Graph
			The graph for which the clustering shall be generated

		Returns:
		--------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeSingletonClustering(G._this))
	def makeOneClustering(self, Graph G):
		"""  Generate a clustering with one cluster consisting of all nodes

		Parameters:
		-----------
		G : networkit.Graph
			The graph for which the clustering shall be generated

		Returns:
		--------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeOneClustering(G._this))
	def makeRandomClustering(self, Graph G, count k):
		"""  Generate a clustering with `k` clusters to which nodes are assigned randomly

		Parameters:
		-----------
		G : networkit.Graph
			The graph for which the clustering shall be generated
		k: count
			The number of clusters that shall be generated

		Returns:
		--------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeRandomClustering(G._this, k))
	def makeContinuousBalancedClustering(self, Graph G, count k):
		"""  Generate a clustering with `k` clusters to which nodes are assigned in continuous blocks

		Parameters:
		-----------
		G : networkit.Graph
			The graph for which the clustering shall be generated
		k: count
			The number of clusters that shall be generated

		Returns:
		--------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeContinuousBalancedClustering(G._this, k))
	def makeNoncontinuousBalancedClustering(self, Graph G, count k):
		"""  Generate a clustering with `k` clusters, the ith node is assigned to cluster i % k. This means that
		for k**2 nodes, this clustering is complementary to the continuous clustering in the sense that no pair
		of nodes that is in the same cluster in one of the clusterings is in the same cluster in the other clustering.

		Parameters:
		-----------
		G : networkit.Graph
			The graph for which the clustering shall be generated
		k: count
			The number of clusters that shall be generated

		Returns:
		--------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeNoncontinuousBalancedClustering(G._this, k))

cdef extern from "<networkit/community/GraphClusteringTools.hpp>" namespace "NetworKit::GraphClusteringTools":

	float getImbalance(_Partition zeta) except +
	_Graph communicationGraph(_Graph graph, _Partition zeta) except +
	count weightedDegreeWithCluster(_Graph graph, _Partition zeta, node u, index cid)
	bool_t isProperClustering(_Graph G, _Partition zeta)
	bool_t isSingletonClustering(_Graph G, _Partition zeta)
	bool_t isOneClustering(_Graph G, _Partition zeta)
	bool_t equalClusterings(_Partition zeta, _Partition eta, _Graph G)

cdef class GraphClusteringTools:
	@staticmethod
	def getImbalance(Partition zeta):
		return getImbalance(zeta._this)
	@staticmethod
	def communicationGraph(Graph graph, Partition zeta):
		return Graph().setThis(communicationGraph(graph._this, zeta._this))
	@staticmethod
	def weightedDegreeWithCluster(Graph graph, Partition zeta, node u, index cid):
		return weightedDegreeWithCluster(graph._this, zeta._this, u, cid)
	@staticmethod
	def isProperClustering(Graph G, Partition zeta):
		return isProperClustering(G._this, zeta._this)
	@staticmethod
	def isSingletonClustering(Graph G, Partition zeta):
		return isSingletonClustering(G._this, zeta._this)
	@staticmethod
	def isOneClustering(Graph G, Partition zeta):
		return isOneClustering(G._this, zeta._this)
	@staticmethod
	def equalClustering(Partition zeta, Partition eta, Graph G):
		return equalClusterings(zeta._this, eta._this, G._this)

cdef extern from "<networkit/community/PartitionIntersection.hpp>":

	cdef cppclass _PartitionIntersection "NetworKit::PartitionIntersection":
		_PartitionIntersection() except +
		_Partition calculate(_Partition zeta, _Partition eta) except +

cdef class PartitionIntersection:
	""" Class for calculating the intersection of two partitions, i.e. the clustering with the fewest clusters
	such that each cluster is a subset of a cluster in both partitions.
	"""
	cdef _PartitionIntersection _this
	def calculate(self, Partition zeta, Partition eta):
		"""  Calculate the intersection of two partitions `zeta` and `eta`

		Parameters:
		-----------
		zeta: networkit.Partition
			The first partition
		eta: networkit.Partition
			The second partition

		Returns:
		--------
		networkit.Partition
			The intersection of zeta and eta
		"""
		return Partition().setThis(self._this.calculate(zeta._this, eta._this))

cdef extern from "<networkit/community/Coverage.hpp>":

	cdef cppclass _Coverage "NetworKit::Coverage":
		_Coverage() except +
		double getQuality(_Partition _zeta, _Graph _G) except +

cdef class Coverage:
	""" Coverage is the fraction of intra-community edges """
	cdef _Coverage _this

	def getQuality(self, Partition zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)


cdef extern from "<networkit/community/EdgeCut.hpp>":

	cdef cppclass _EdgeCut "NetworKit::EdgeCut":
		_EdgeCut() except +
		double getQuality(_Partition _zeta, _Graph _G) except +

cdef class EdgeCut:
	""" Edge cut is the total weight of inter-community edges"""
	cdef _EdgeCut _this

	def getQuality(self, Partition zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)


cdef extern from "<networkit/community/Modularity.hpp>":

	cdef cppclass _Modularity "NetworKit::Modularity":
		_Modularity() except +
		double getQuality(_Partition _zeta, _Graph _G) nogil except +


cdef class Modularity:
	"""	Modularity is a quality index for community detection.
	It assigns a quality value in [-0.5, 1.0] to a partition of a graph which is higher for more modular networks and
	partitions which better capture the modular structure. See also http://en.wikipedia.org/wiki/Modularity_(networks).

 	Notes
	-----
	Modularity is defined as:

	.. math:: mod(\zeta) := \\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)} - \\frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }

	"""
	cdef _Modularity _this

	def getQuality(self, Partition zeta, Graph G):
		cdef double ret
		with nogil:
			ret = self._this.getQuality(zeta._this, G._this)
		return ret

cdef extern from "<networkit/community/HubDominance.hpp>":

	cdef cppclass _HubDominance "NetworKit::HubDominance":
		_HubDominance() except +
		double getQuality(_Partition _zeta, _Graph _G) except +
		double getQuality(_Cover _zeta, _Graph _G) except +

cdef class HubDominance:
	"""
	A quality measure that measures the dominance of hubs in clusters. The hub dominance of a single
	cluster is defined as the maximum cluster-internal degree of a node in that cluster divided by
	the maximum cluster-internal degree, i.e. the number of nodes in the cluster minus one. The
	value for all clusters is defined as the average of all clusters.

	Strictly speaking this is not a quality measure as this is rather dependent on the type of the
	considered graph, for more information see
	Lancichinetti A, Kivel M, Saramki J, Fortunato S (2010)
	Characterizing the Community Structure of Complex Networks
	PLoS ONE 5(8): e11976. doi: 10.1371/journal.pone.0011976
	http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011976
	"""

	cdef _HubDominance _this

	def getQuality(self, PartitionCover zeta, Graph G):
		"""
		Calculates the dominance of hubs in the given Partition or Cover of the given
		Graph.

		Parameters:
		-----------
		zeta : networkit.Partition or networkit.Cover
			The Partition or Cover for which the hub dominance shall be calculated
		G : networkit.Graph
			The Graph to which zeta belongs

		Returns:
		--------
		double
			The average hub dominance in the given Partition or Cover
		"""
		return self._this.getQuality(zeta._this, G._this)

cdef extern from "<networkit/community/PLM.hpp>":

	cdef cppclass _PLM "NetworKit::PLM"(_CommunityDetectionAlgorithm):
		_PLM(_Graph _G) except +
		_PLM(_Graph _G, bool_t refine, double gamma, string par, count maxIter, bool_t turbo, bool_t recurse, _Partition _zeta) except +
		map[string, vector[count]] getTiming() except +

cdef extern from "<networkit/community/PLM.hpp>" namespace "NetworKit::PLM":

	pair[_Graph, vector[node]] PLM_coarsen "NetworKit::PLM::coarsen" (const _Graph& G, const _Partition& zeta) except +
	_Partition PLM_prolong "NetworKit::PLM::prolong"(const _Graph& Gcoarse, const _Partition& zetaCoarse, const _Graph& Gfine, vector[node] nodeToMetaNode) except +


cdef class PLM(CommunityDetector):
	""" Parallel Louvain Method - the Louvain method, optionally extended to
		a full multi-level algorithm with refinement

		Parameters:
		-----------
		G : networkit.Graph
			A graph.
		refine : bool, optional
			Add a second move phase to refine the communities.
		gamma : double
			Multi-resolution modularity parameter:
			1.0 -> standard modularity
	 		0.0 -> one community
	 		2m 	-> singleton communities
		par : string
			parallelization strategy
		maxIter : count
			maximum number of iterations for move phase
		turbo : bool, optional
			faster but uses O(n) additional memory per thread
		recurse: bool, optional
			use recursive coarsening, see http://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.049902 for some explanations (default: true)
	"""

	def __cinit__(self, Graph G not None, refine=False, gamma=1.0, par="balanced", maxIter=32, turbo=True, recurse=True, Partition zeta=Partition(0)):
		self._G = G
		self._this = new _PLM(G._this, refine, gamma, stdstring(par), maxIter, turbo, recurse, zeta._this)

	def getTiming(self):
		"""  Get detailed time measurements.
		"""
		return (<_PLM*>(self._this)).getTiming()

	@staticmethod
	def coarsen(Graph G, Partition zeta, bool_t parallel = False):
		cdef pair[_Graph, vector[node]] result = move(PLM_coarsen(G._this, zeta._this))
		return (Graph().setThis(result.first), result.second)

	@staticmethod
	def prolong(Graph Gcoarse, Partition zetaCoarse, Graph Gfine, vector[node] nodeToMetaNode):
		return Partition().setThis(PLM_prolong(Gcoarse._this, zetaCoarse._this, Gfine._this, nodeToMetaNode))

cdef extern from "<networkit/community/LouvainMapEquation.hpp>":
	cdef cppclass _LouvainMapEquation "NetworKit::LouvainMapEquation"(_CommunityDetectionAlgorithm):
		_LouvainMapEquation(_Graph, bool, count, string ) except +

cdef class LouvainMapEquation(CommunityDetector):
	"""
	Community detection algorithm based on the Louvain algorithm. Uses the Map Equation to find communities.

	Parameters
	----------
	G : networkit.Graph
		The graph on which the algorithm has to run.
	hierarchical: bool
		(optional) Iteratively create a graph of the locally optimal clusters and optimize locally on that graph.
	maxIterations: count
		(optional) The maximum number of local move iterations.
	parallelizationStrategy: string
		(optional) relaxmap, synchronous, or none. default relaxmap.
	"""

	def __cinit__(self, Graph G not None, hierarchical = False, maxIterations = 32, parallelizationStrategy = "relaxmap"):
		self._G = G
		self._this = new _LouvainMapEquation(G._this, hierarchical, maxIterations, stdstring(parallelizationStrategy))

cdef extern from "<networkit/community/PLP.hpp>":

	cdef cppclass _PLP "NetworKit::PLP"(_CommunityDetectionAlgorithm):
		_PLP(_Graph _G, count updateThreshold, count maxIterations) except +
		_PLP(_Graph _G, _Partition baseClustering, count updateThreshold) except +
		count numberOfIterations() except +
		vector[count] getTiming() except +


cdef class PLP(CommunityDetector):
	""" Parallel label propagation for community detection:
	Moderate solution quality, very short time to solution.

	Parameters:
	-----------
	G : networkit.Graph
		The graph on which the algorithm has to run.
	updateThreshold : int
		number of nodes that have to be changed in each iteration so that a new iteration starts.
	baseClustering : networkit.Partition
		PLP needs a base clustering to start from; if none is given the algorithm will run on a singleton clustering.

	Notes
	-----
	As described in Ovelgoenne et al: An Ensemble Learning Strategy for Graph Clustering
 	Raghavan et al. proposed a label propagation algorithm for graph clustering.
 	This algorithm initializes every vertex of a graph with a unique label. Then, in iterative
 	sweeps over the set of vertices the vertex labels are updated. A vertex gets the label
 	that the maximum number of its neighbors have. The procedure is stopped when every vertex
 	has the label that at least half of its neighbors have.
	"""

	def __cinit__(self, Graph G not None, count updateThreshold=none, count maxIterations=none, Partition baseClustering=None,):
		"""
		Constructor to the Parallel label propagation community detection algorithm.

		"""
		self._G = G


		if baseClustering is None:
			self._this = new _PLP(G._this, updateThreshold, maxIterations)
		else:
			self._this = new _PLP(G._this, baseClustering._this, updateThreshold)


	def numberOfIterations(self):
		""" Get number of iterations in last run.

		Returns:
		--------
		count
			The number of iterations.
		"""
		return (<_PLP*>(self._this)).numberOfIterations()

	def getTiming(self):
		""" Get list of running times for each iteration.

		Returns:
		--------
		count
			The list of running times in milliseconds.
		"""
		return (<_PLP*>(self._this)).getTiming()

cdef extern from "<networkit/community/LPDegreeOrdered.hpp>":

	cdef cppclass _LPDegreeOrdered "NetworKit::LPDegreeOrdered"(_CommunityDetectionAlgorithm):
		_LPDegreeOrdered(_Graph _G) except +
		count numberOfIterations()

cdef class LPDegreeOrdered(CommunityDetector):
	""" Label propagation-based community detection algorithm which processes nodes in increasing order of node degree.	"""

	def __cinit__(self, Graph G not None):
		self._G = G
		self._this = new _LPDegreeOrdered(G._this)

	def numberOfIterations(self):
		""" Get number of iterations in last run.

		Returns:
		--------
		count
			Number of iterations.
		"""
		return (<_LPDegreeOrdered*>(self._this)).numberOfIterations()

cdef extern from "<networkit/community/CutClustering.hpp>":

	cdef cppclass _CutClustering "NetworKit::CutClustering"(_CommunityDetectionAlgorithm):
		_CutClustering(_Graph _G) except +
		_CutClustering(_Graph _G, edgeweight alpha) except +

cdef extern from "<networkit/community/CutClustering.hpp>" namespace "NetworKit::CutClustering":

	map[double, _Partition] CutClustering_getClusterHierarchy "NetworKit::CutClustering::getClusterHierarchy"(const _Graph& G) nogil except +


cdef class CutClustering(CommunityDetector):
	"""
	Cut clustering algorithm as defined in
	Flake, Gary William; Tarjan, Robert E.; Tsioutsiouliklis, Kostas. Graph Clustering and Minimum Cut Trees.
	Internet Mathematics 1 (2003), no. 4, 385--408.

	Parameters:
	-----------
	G : networkit.Graph
	alpha : double
		The parameter for the cut clustering algorithm
	"""
	def __cinit__(self, Graph G not None,  edgeweight alpha):
		self._G = G
		self._this = new _CutClustering(G._this, alpha)

	@staticmethod
	def getClusterHierarchy(Graph G not None):
		""" Get the complete hierarchy with all possible parameter values.

		Each reported parameter value is the lower bound for the range in which the corresponding clustering is calculated by the cut clustering algorithm.

		Warning: all reported parameter values are slightly too high in order to avoid wrong clusterings because of numerical inaccuracies.
		Furthermore the completeness of the hierarchy cannot be guaranteed because of these inaccuracies.
		This implementation hasn't been optimized for performance.

		Parameters:
		-----------
		G : networkit.Graph
			The graph.

		Returns:
		--------
		dict
			A dictionary with the parameter values as keys and the corresponding Partition instances as values
		"""
		cdef map[double, _Partition] result
		# FIXME: this probably copies the whole hierarchy because of exception handling, using move might fix this
		with nogil:
			result = CutClustering_getClusterHierarchy(G._this)
		pyResult = {}
		# FIXME: this code copies the partitions a lot!
		for res in result:
			pyResult[res.first] = Partition().setThis(res.second)
		return pyResult

cdef class DissimilarityMeasure:
	""" Abstract base class for partition/community dissimilarity measures """
	# TODO: use conventional class design of parametrized constructor, run-method and getters
	pass

cdef extern from "<networkit/community/NodeStructuralRandMeasure.hpp>":

	cdef cppclass _NodeStructuralRandMeasure "NetworKit::NodeStructuralRandMeasure":
		_NodeStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class NodeStructuralRandMeasure(DissimilarityMeasure):
	""" The node-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering all pairs of nodes.
	"""
	cdef _NodeStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret


cdef extern from "<networkit/community/GraphStructuralRandMeasure.hpp>":

	cdef cppclass _GraphStructuralRandMeasure "NetworKit::GraphStructuralRandMeasure":
		_GraphStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class GraphStructuralRandMeasure(DissimilarityMeasure):
	""" The graph-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering connected pairs of nodes.
	"""
	cdef _GraphStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret


cdef extern from "<networkit/community/JaccardMeasure.hpp>":

	cdef cppclass _JaccardMeasure "NetworKit::JaccardMeasure":
		_JaccardMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class JaccardMeasure(DissimilarityMeasure):
	""" TODO:
	"""
	cdef _JaccardMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret

cdef extern from "<networkit/community/NMIDistance.hpp>":

	cdef cppclass _NMIDistance "NetworKit::NMIDistance":
		_NMIDistance() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class NMIDistance(DissimilarityMeasure):
	""" The NMI distance assigns a similarity value in [0,1] to two partitions
		of a graph.
	"""
	cdef _NMIDistance _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret

cdef extern from "<networkit/community/AdjustedRandMeasure.hpp>":

	cdef cppclass _AdjustedRandMeasure "NetworKit::AdjustedRandMeasure":
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class AdjustedRandMeasure(DissimilarityMeasure):
	"""
	The adjusted rand dissimilarity measure as proposed by Huber and Arabie in "Comparing partitions" (http://link.springer.com/article/10.1007/BF01908075)
	"""
	cdef _AdjustedRandMeasure _this

	def getDissimilarity(self, Graph G not None, Partition first not None, Partition second not None):
		"""
		Get the adjust rand dissimilarity. Runs in O(n log(n)).

		Note that the dissimilarity can be larger than 1 if the partitions are more different than expected in the random model.

		Parameters:
		-----------
		G : networkit.Graph
			The graph on which the partitions shall be compared
		zeta : networkit.Partition
			The first partiton
		eta : networkit.Partition
			The second partition

		Returns:
		--------
		double
			The adjusted rand dissimilarity
		"""
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret

cdef extern from "<networkit/community/LocalCommunityEvaluation.hpp>":

	cdef cppclass _LocalCommunityEvaluation "NetworKit::LocalCommunityEvaluation"(_Algorithm):
		double getWeightedAverage() except +
		double getUnweightedAverage() except +
		double getMaximumValue() except +
		double getMinimumValue() except +
		double getValue(index i) except +
		vector[double] getValues() except +
		bool_t isSmallBetter() except +

cdef class LocalCommunityEvaluation(Algorithm):
	"""
	Virtual base class of all evaluation methods for a single clustering which is based on the evaluation of single clusters.
	This is the base class both for Partitions as well as for Covers.
	"""
	def __init__(self, *args, **namedargs):
		if type(self) == LocalCommunityEvaluation:
			raise RuntimeError("Error, you may not use LocalCommunityEvaluation directly, use a sub-class instead")

	def getWeightedAverage(self):
		""" Get the average value weighted by cluster size.

		Returns:
		--------
		double:
			The weighted average value.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getWeightedAverage()

	def getUnweightedAverage(self):
		""" Get the (unweighted) average value of all clusters.

		Returns:
		--------
		double:
			The unweighted average value.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getUnweightedAverage()

	def getMaximumValue(self):
		""" Get the maximum value of all clusters.

		Returns:
		--------
		double:
			The maximum value.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getMaximumValue()

	def getMinimumValue(self):
		""" Get the minimum value of all clusters.

		Returns:
		--------
		double:
			The minimum value.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getMinimumValue()

	def getValue(self, index i):
		""" Get the value of the specified cluster.

		Parameters:
		-----------
		i : index
			The cluster to get the value for.

		Returns:
		--------
		double:
			The value of cluster i.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getValue(i)

	def getValues(self):
		""" Get the values of all clusters.

		Returns:
		--------
		list[double]:
			The values of all clusters.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getValues()

	def isSmallBetter(self):
		""" If small values are better (otherwise large values are better).

		Returns:
		--------
		bool:
			If small values are better.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).isSmallBetter()

cdef extern from "<networkit/community/LocalPartitionEvaluation.hpp>":

	cdef cppclass _LocalPartitionEvaluation "NetworKit::LocalPartitionEvaluation"(_LocalCommunityEvaluation):
		pass

cdef class LocalPartitionEvaluation(LocalCommunityEvaluation):
	"""
	Virtual base class of all evaluation methods for a single clustering which is based on the evaluation of single clusters.
	This is the base class for Partitions.
	"""
	cdef Graph _G
	cdef Partition _P

	def __init__(self, *args, **namedargs):
		if type(self) == LocalPartitionEvaluation:
			raise RuntimeError("Error, you may not use LocalPartitionEvaluation directly, use a sub-class instead")

	def __cinit__(self, Graph G not None, Partition P not None, *args, **namedargs):
		self._G = G
		self._P = P

	def __dealloc__(self):
		# Just to be sure that everything is properly deleted
		self._G = None
		self._P = None


cdef extern from "<networkit/community/LocalCoverEvaluation.hpp>":

	cdef cppclass _LocalCoverEvaluation "NetworKit::LocalCoverEvaluation"(_LocalCommunityEvaluation):
		pass


cdef class LocalCoverEvaluation(LocalCommunityEvaluation):
	"""
	Virtual base class of all evaluation methods for a single clustering which is based on the evaluation of single clusters.
	This is the base class for Covers.
	"""
	cdef Graph _G
	cdef Cover _C

	def __init__(self, *args, **namedargs):
		if type(self) == LocalCoverEvaluation:
			raise RuntimeError("Error, you may not use LocalCoverEvaluation directly, use a sub-class instead")

	def __cinit__(self, Graph G not None, Cover C not None, *args, **namedargs):
		self._G = G
		self._C = C

	def __dealloc__(self):
		# Just to be sure that everything is properly deleted
		self._G = None
		self._C = None

cdef extern from "<networkit/community/IntrapartitionDensity.hpp>":

	cdef cppclass _IntrapartitionDensity "NetworKit::IntrapartitionDensity"(_LocalPartitionEvaluation):
		_IntrapartitionDensity(_Graph G, _Partition P) except +
		double getGlobal() except +

cdef class IntrapartitionDensity(LocalPartitionEvaluation):
	"""
	The intra-cluster density of a partition is defined as the number of existing edges divided by the number of possible edges.
	The global value is the sum of all existing intra-cluster edges divided by the sum of all possible intra-cluster edges.

	Parameters:
	-----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _IntrapartitionDensity(self._G._this, self._P._this)

	def getGlobal(self):
		""" Get the global intra-cluster density.

		Returns:
		--------
		double:
			The global intra-cluster density.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_IntrapartitionDensity*>(self._this)).getGlobal()


cdef extern from "<networkit/community/IsolatedInterpartitionConductance.hpp>":

	cdef cppclass _IsolatedInterpartitionConductance "NetworKit::IsolatedInterpartitionConductance"(_LocalPartitionEvaluation):
		_IsolatedInterpartitionConductance(_Graph G, _Partition P) except +

cdef class IsolatedInterpartitionConductance(LocalPartitionEvaluation):
	"""
	Isolated inter-partition conductance is a measure for how well a partition
	(communtiy/cluster) is separated from the rest of the graph.

	The conductance of a partition is defined as the weight of the cut divided
	by the volume (the sum of the degrees) of the nodes in the partition or the
	nodes in the rest of the graph, whatever is smaller. Small values thus indicate
	that the cut is small compared to the volume of the smaller of the separated
	parts. For the whole partitions usually the maximum or the unweighted average
	is used.

	See also Experiments on Density-Constrained Graph Clustering,
	Robert Grke, Andrea Kappes and  Dorothea Wagner, JEA 2015:
	http://dx.doi.org/10.1145/2638551

	Parameters:
	-----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _IsolatedInterpartitionConductance(self._G._this, self._P._this)

cdef extern from "<networkit/community/IsolatedInterpartitionExpansion.hpp>":

	cdef cppclass _IsolatedInterpartitionExpansion "NetworKit::IsolatedInterpartitionExpansion"(_LocalPartitionEvaluation):
		_IsolatedInterpartitionExpansion(_Graph G, _Partition P) except +

cdef class IsolatedInterpartitionExpansion(LocalPartitionEvaluation):
	"""
	Isolated inter-partition expansion is a measure for how well a partition
	(communtiy/cluster) is separated from the rest of the graph.

	The expansion of a partition is defined as the weight of the cut divided
	by number of nodes in the partition or in the rest of the graph, whatever
	is smaller. Small values thus indicate that the cut is small compared to
	the size of the smaller of the separated parts. For the whole partitions
	usually the maximum or the unweighted average is used. Note that expansion
	values can be larger than 1.

	See also Experiments on Density-Constrained Graph Clustering,
	Robert Grke, Andrea Kappes and Dorothea Wagner, JEA 2015:
	http://dx.doi.org/10.1145/2638551

	Parameters:
	-----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _IsolatedInterpartitionExpansion(self._G._this, self._P._this)

cdef extern from "<networkit/community/CoverHubDominance.hpp>":

	cdef cppclass _CoverHubDominance "NetworKit::CoverHubDominance"(_LocalCoverEvaluation):
		_CoverHubDominance(_Graph G, _Cover C) except +

cdef class CoverHubDominance(LocalCoverEvaluation):
	"""
	A quality measure that measures the dominance of hubs in clusters. The hub dominance of a single
	cluster is defined as the maximum cluster-internal degree of a node in that cluster divided by
	the maximum cluster-internal degree, i.e. the number of nodes in the cluster minus one. The
	value for all clusters is defined as the average of all clusters.
	This implementation is a natural generalization of this measure for covers.
	Strictly speaking this is not a quality measure as this is rather dependent on the type of the
	considered graph, for more information see
	Lancichinetti A, Kivel M, Saramki J, Fortunato S (2010)
	Characterizing the Community Structure of Complex Networks
	PLoS ONE 5(8): e11976. doi: 10.1371/journal.pone.0011976
	http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011976

	Parameters:
	-----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	C : networkit.Cover
		The cover that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _CoverHubDominance(self._G._this, self._C._this)

cdef extern from "<networkit/community/PartitionHubDominance.hpp>":

	cdef cppclass _PartitionHubDominance "NetworKit::PartitionHubDominance"(_LocalPartitionEvaluation):
		_PartitionHubDominance(_Graph G, _Partition C) except +

cdef class PartitionHubDominance(LocalPartitionEvaluation):
	"""
	A quality measure that measures the dominance of hubs in clusters. The hub dominance of a single
	cluster is defined as the maximum cluster-internal degree of a node in that cluster divided by
	the maximum cluster-internal degree, i.e. the number of nodes in the cluster minus one. The
	value for all clusters is defined as the average of all clusters.
	Strictly speaking this is not a quality measure as this is rather dependent on the type of the
	considered graph, for more information see
	Lancichinetti A, Kivel M, Saramki J, Fortunato S (2010)
	Characterizing the Community Structure of Complex Networks
	PLoS ONE 5(8): e11976. doi: 10.1371/journal.pone.0011976
	http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011976

	Parameters:
	-----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _PartitionHubDominance(self._G._this, self._P._this)

cdef extern from "<networkit/community/PartitionFragmentation.hpp>":

	cdef cppclass _PartitionFragmentation "NetworKit::PartitionFragmentation"(_LocalPartitionEvaluation):
		_PartitionFragmentation(_Graph G, _Partition C) except +

cdef class PartitionFragmentation(LocalPartitionEvaluation):
	"""
	This measure evaluates how fragmented a partition is. The fragmentation of a single cluster is defined as one minus the
	number of nodes in its maximum connected componented divided by its total number of nodes. Smaller values thus indicate a smaller fragmentation.

	Parameters:
	-----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _PartitionFragmentation(self._G._this, self._P._this)

cdef extern from "<networkit/community/StablePartitionNodes.hpp>":

	cdef cppclass _StablePartitionNodes "NetworKit::StablePartitionNodes"(_LocalPartitionEvaluation):
		_StablePartitionNodes(_Graph G, _Partition C) except +
		bool_t isStable(node u) except +

cdef class StablePartitionNodes(LocalPartitionEvaluation):
	"""
	Evaluates how stable a given partition is. A node is considered to be stable if it has strictly more connections
	to its own partition than to other partitions. Isolated nodes are considered to be stable.
	The value of a cluster is the percentage of stable nodes in the cluster.
	Larger values indicate that a clustering is more stable and thus better defined.

	Parameters:
	-----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _StablePartitionNodes(self._G._this, self._P._this)


	def isStable(self, node u):
		"""
		Check if a given node is stable, i.e. more connected to its own partition than to other partitions.

		Parameters:
		-----------
		u : node
			The node to check

		Returns:
		--------
		bool
			If the node u is stable.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_StablePartitionNodes*>(self._this)).isStable(u)


cdef extern from "<networkit/community/CoverF1Similarity.hpp>":

	cdef cppclass _CoverF1Similarity "NetworKit::CoverF1Similarity"(_LocalCoverEvaluation):
		_CoverF1Similarity(_Graph G, _Cover C, _Cover reference) except +

cdef class CoverF1Similarity(LocalCoverEvaluation):
	"""
	Compare a given cover to a reference cover using the F1 measure.
	This is a typical similarity measure used to compare the found
	overlapping community structure to a ground truth community
	structure. Each cluster is compared to the best-matching reference
	cluster (in terms of highest F1 score). A value of 1 indicates
	perfect agreement while a while of 0 indicates complete
	disagreement. An example where this measure is used is the
	following paper:

	Alessandro Epasto, Silvio Lattanzi, and Renato Paes
	Leme. 2017. Ego-Splitting Framework: from Non-Overlapping to
	Overlapping Clusters. In Proceedings of the 23rd ACM SIGKDD
	International Conference on Knowledge Discovery and Data Mining
	(KDD '17). ACM, New York, NY, USA, 145-154. DOI:
	https://doi.org/10.1145/3097983.3098054

	Parameters:
	-----------
	G : Graph
		The graph on which the evaluation is performed.
	C : Cover
		The cover that shall be evaluated
        reference : Cover
		The cover to which the similarity shall be computed
	"""
	cdef Cover _reference
	def __cinit__(self, Graph G not None, Cover C not None, Cover reference not None):
		self._this = new _CoverF1Similarity(G._this, C._this, reference._this)
		self._reference = reference
		assert(self._G == G)
		assert(self._C == C)

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
	print("Communities detected in {:.5f} [s]".format(t.elapsed))
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
		raise RuntimeError("There is no core for the specified k")

	C = graph.Subgraph().fromNodes(G, kCore)	# FIXME: node indices are not preserved

	#properties.overview(C)

	return detectCommunities(C, algo, inspect)
"""
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
"""


cdef extern from "<networkit/community/OverlappingNMIDistance.hpp>" namespace "NetworKit::OverlappingNMIDistance":

	cdef enum _Normalization "NetworKit::OverlappingNMIDistance::Normalization":
		MIN,
		GEOMETRIC_MEAN,
		ARITHMETIC_MEAN,
		MAX,
		JOINT_ENTROPY

cdef extern from "<networkit/community/OverlappingNMIDistance.hpp>":

	cdef cppclass _OverlappingNMIDistance "NetworKit::OverlappingNMIDistance":
		_OverlappingNMIDistance() except +
		_OverlappingNMIDistance(_Normalization normalization) except +
		void setNormalization(_Normalization normalization) except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +
		double getDissimilarity(_Graph G, _Cover first, _Cover second) nogil except +

cdef class OverlappingNMIDistance(DissimilarityMeasure):
	"""
	Compare two covers using the overlapping normalized mutual information measure. This is a dissimilarity measure with
	a range of [0, 1]. A value of 0 indicates a perfect agreement while a 1 indicates complete disagreement.

	For the `OverlappingNMIDistance.Max` normalization, this is the measure introduced in [NMI13]. Other normalization
	methods result in similar measures.

	Parameters
	----------
	normalization : {Min, GeometricMean, ArithmeticMean, Max, JointEntropy}, optional
		The default is OverlappingNMIDistance.Max.

	Raises
	------
	ValueError
	    If `normalization` is not one of the available methods.

	References
	----------
	[NMI13]
		McDaid, Aaron F., Derek Greene, and Neil Hurley. "Normalized Mutual Information to Evaluate Overlapping
		Community Finding Algorithms." ArXiv:1110.2515 [Physics], August 2, 2013. http://arxiv.org/abs/1110.2515.
	"""
	cdef _OverlappingNMIDistance _this

	Min = _Normalization.MIN
	GeometricMean = _Normalization.GEOMETRIC_MEAN
	ArithmeticMean = _Normalization.ARITHMETIC_MEAN
	Max = _Normalization.MAX
	JointEntropy = _Normalization.JOINT_ENTROPY

	def __cinit__(self, _Normalization normalization = _Normalization.MAX):
		self._validateNormalization(normalization)
		self._this = _OverlappingNMIDistance(normalization)

	def setNormalization(self, _Normalization normalization):
		"""
		Set the normalization method.

		Parameters
		----------
		normalization : {Min, GeometricMean, ArithmeticMean, Max, JointEntropy}

		Raises
		------
		ValueError
		    If `normalization` is not one of the available methods.
		"""
		self._validateNormalization(normalization)
		self._this.setNormalization(normalization)

	def getDissimilarity(self, Graph G, PartitionCover first, PartitionCover second):
		"""
		Calculate the dissimilarity.

		Parameters
		----------
		G : networkit.Graph
		first : networkit.Partition or networkit.Cover
		second : networkit.Partition or networkit.Cover
			Must be the same type as `first`.

		Raises
		------
		TypeError
		    If `first` and `second` do not have the same type.
		ValueError
		    If `G`, `first` and `second` do not have the matching number of nodes.

		Returns
		-------
		distance : float
		"""
		cdef double ret
		if isinstance(first, Partition) and isinstance(second, Partition):
			with nogil:
				ret = self._this.getDissimilarity(G._this, (<Partition>(first))._this, (<Partition>(second))._this)
		elif isinstance(first, Cover) and isinstance(second, Cover):
			with nogil:
				ret = self._this.getDissimilarity(G._this, (<Cover>(first))._this, (<Cover>(second))._this)
		else:
			raise TypeError("Error, first and second must both be either a Partition or a Cover")
		return ret

	def _validateNormalization(self, _Normalization normalization):
		if normalization not in {OverlappingNMIDistance.Min, OverlappingNMIDistance.GeometricMean,
				OverlappingNMIDistance.ArithmeticMean, OverlappingNMIDistance.Max, OverlappingNMIDistance.JointEntropy}:
			raise ValueError("Error, invalid normalization method")
