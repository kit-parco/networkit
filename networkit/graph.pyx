# distutils: language=c++

from cython.operator import dereference, preincrement

from .base import Algorithm
from .helpers import stdstring, pystring
from .traversal import Traversal
from .graphio import readGraph, writeGraph, Format
from fs.memoryfs import MemoryFS
import os

cdef class Graph:

	""" An undirected graph (with optional weights) and parallel iterator methods.

		Graph(n=0, weighted=False, directed=False)

		Create a graph of `n` nodes. The graph has assignable edge weights if `weighted` is set to True.
		If `weighted` is set to False each edge has edge weight 1.0 and any other weight assignment will
		be ignored.

	    Parameters:
	    -----------
		n : count, optional
			Number of nodes.
		weighted : bool, optional
			If set to True, the graph can have edge weights other than 1.0.
		directed : bool, optional
			If set to True, the graph will be directed.
	"""

	def __cinit__(self, n=0, bool_t weighted=False, bool_t directed=False):
		if isinstance(n, Graph):
			self._this = move(_Graph((<Graph>n)._this, weighted, directed))
		else:
			self._this = move(_Graph(<count>n, weighted, directed))

	cdef setThis(self, _Graph& other):
		swap[_Graph](self._this, other)
		return self

	def __copy__(self):
		"""
		Generates a copy of the graph
		"""
		return Graph().setThis(_Graph(self._this))

	def __deepcopy__(self, memo):
		"""
		Generates a (deep) copy of the graph
		"""
		return Graph().setThis(_Graph(self._this))

	def __str__(self):
		return "NetworKit.Graph(n={0}, m={1})".format(self.numberOfNodes(), self.numberOfEdges())
	
	def __getstate__(self):
		if not os.path.isdir('./tmp/'):
			os.makedirs('./tmp')
		i = 0
		filename = "./tmp/tmpGraph" + str(i)
		while os.path.exists(filename):
			i = i+1
			filename = "./tmp/tmpGraph" + str(i)
		writeGraph(self, filename, Format.NetworkitBinary, chunks = 32, weightsType = 5, preserveEdgeIndex = True)
		
		file = open(filename, mode='rb')
		graph_bytes = file.read()
		file.close()
		os.remove(filename)
		
		return graph_bytes
	
	def __setstate__(self, state):
		if not os.path.isdir('./tmp/'):
			os.makedirs('./tmp')
		i = 0
		filename = "./tmp/tmpGraph" + str(i)
		while os.path.exists(filename):
			i = i+1
			filename = "./tmp/tmpGraph" + str(i)
		
		file = open(filename, mode='wb')
		file.write(state)
		file.close()
		os.remove(filename)
	
		self.this = readGraph(state, Format.NetworkitBinary)

	def indexEdges(self, bool_t force = False):
		"""
		Assign integer ids to edges.

		Parameters:
		-----------
		force : bool
			Force re-indexing of edges.

		"""
		self._this.indexEdges(force)

	def hasEdgeIds(self):
		"""
		Returns true if edges have been indexed

		Returns:
		--------
		bool
			If edges have been indexed
		"""
		return self._this.hasEdgeIds()

	def edgeId(self, node u, node v):
		"""
		Returns:
		--------
		edgeid
			id of the edge
		"""
		return self._this.edgeId(u, v)

	def numberOfNodes(self):
		"""
		Get the number of nodes in the graph.

		Returns:
		--------
		count
			The number of nodes.
		"""
		return self._this.numberOfNodes()

	def numberOfEdges(self):
		"""
		Get the number of edges in the graph.

		Returns:
	 	--------
		count
			The number of edges.
		"""
		return self._this.numberOfEdges()

	def upperNodeIdBound(self):
		"""
		Get an upper bound for the node ids in the graph

		Returns:
		--------
		count
			An upper bound for the node ids in the graph
		"""
		return self._this.upperNodeIdBound()

	def upperEdgeIdBound(self):
		"""
		Get an upper bound for the edge ids in the graph

		Returns:
		--------
		count
			An upper bound for the edge ids in the graph
		"""
		return self._this.upperEdgeIdBound()

	def degree(self, u):
		"""
		Get the number of neighbors of `v`.

		Parameters:
		-----------
		v : node
			Node.

		Returns:
		--------
		count
			The number of neighbors.
		"""
		return self._this.degree(u)

	def degreeIn(self, u):
		return self._this.degreeIn(u)

	def degreeOut(self, u):
		return self._this.degreeOut(u)

	def weightedDegree(self, u, countSelfLoopsTwice=False):
		"""
		Returns the weighted out-degree of u.

		For directed graphs this is the sum of weights of all outgoing edges of u.

		Parameters:
		-----------
		u : node
			Node.
		countSelfLoopsTwice : bool
			If set to true, self-loops will be counted twice

		Returns:
		--------
		double
			The weighted out-degree of u.
		"""
		return self._this.weightedDegree(u, countSelfLoopsTwice)

	def weightedDegreeIn(self, u, countSelfLoopsTwice=False):
		"""
		Returns the weighted in-degree of u.

		For directed graphs this is the sum of weights of all ingoing edges of u.

		Parameters:
		-----------
		u : node
			Node.
		countSelfLoopsTwice : bool
			If set to true, self-loops will be counted twice

		Returns:
		--------
		double
			The weighted in-degree of u.
		"""
		return self._this.weightedDegreeIn(u, countSelfLoopsTwice)

	def isIsolated(self, u):
		"""
		If the node `u` is isolated

		Parameters:
		-----------
		u : node
			Node.

		Returns:
		--------
		bool
			If the node is isolated
		"""
		return self._this.isIsolated(u)

	def addNode(self):
		""" Add a new node to the graph and return it.

		Returns:
		--------
		node
			The new node.
		"""
		return self._this.addNode()

	def addNodes(self, numberOfNewNodes):
		""" Add numberOfNewNodes many new nodes to the graph and return
		the id of the last node added.

		Parameters:
		-----------
		numberOfNewNodes : node
			Number of nodes to be added.

		Returns:
		--------
		node
			The id of the last node added.
		"""
		assert(numberOfNewNodes >= 0)
		return self._this.addNodes(numberOfNewNodes)

	def removeNode(self, u):
		""" Remove a node `v` and all incident edges from the graph.

		Incoming as well as outgoing edges will be removed.

		Parameters:
		-----------
		u : node
			Id of node to be removed.
		"""
		self._this.removeNode(u)

	def restoreNode(self, u):
		""" Restores a previously deleted node `u` with its previous id in the graph.

		Parameters:
		-----------
		u : node
			Node.
		"""
		self._this.restoreNode(u)

	def hasNode(self, u):
		""" Checks if the Graph has the node `u`, i.e. if `u` hasn't been deleted and is in the range of valid ids.

		Parameters:
		-----------
		u : node
			Id of node queried.

		Returns:
		--------
		bool
			If the Graph has the node `u`
		"""
		return self._this.hasNode(u)

	def addEdge(self, u, v, w=1.0, addMissing = False):
		""" Insert an undirected edge between the nodes `u` and `v`. If the graph is weighted you can optionally
		set a weight for this edge. The default weight is 1.0.
		If one or both end-points do not exists and addMissing is set, they are silently added.
		Caution: It is not checked whether this edge already exists, thus it is possible to create multi-edges.

	 	Parameters:
	 	-----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.
		w : edgeweight, optional
			Edge weight.
		addMissing : optional, default: False
			Add missing endpoints if necessary (i.e., increase numberOfNodes).
		"""
		if not (self._this.hasNode(u) and self._this.hasNode(v)):
			if not addMissing:
				raise RuntimeError("Cannot create edge ({0}, {1}) as at least one end point does not exist".format(u,v))

			k = max(u, v)
			if k >= self._this.upperNodeIdBound():
				self._this.addNodes(k - self._this.upperNodeIdBound() + 1)

			if not self._this.hasNode(u):
				self._this.restoreNode(u)

			if not self._this.hasNode(v):
				self._this.restoreNode(v)

		self._this.addEdge(u, v, w)
		return self

	def setWeight(self, u, v, w):
		""" Set the weight of an edge. If the edge does not exist, it will be inserted.

		Parameters:
		-----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.
		w : edgeweight
			Edge weight.
		"""
		self._this.setWeight(u, v, w)
		return self

	def increaseWeight(self, u, v, w):
		""" Increase the weight of an edge. If the edge does not exist, it will be inserted.

		Parameters:
		-----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.
		w : edgeweight
			Edge weight.
		"""
		self._this.increaseWeight(u, v, w)
		return self

	def removeEdge(self, u, v):
		""" Removes the undirected edge {`u`,`v`}.

		Parameters:
		-----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.
		"""
		self._this.removeEdge(u, v)
		return self

	def removeAllEdges(self):
		""" Removes all the edges in the graph.
		"""
		self._this.removeAllEdges()

	def removeSelfLoops(self):
		""" Removes all self-loops from the graph.
		"""
		self._this.removeSelfLoops()

	def removeMultiEdges(self):
		""" Removes all multi-edges from the graph.
		"""
		self._this.removeMultiEdges()

	def swapEdge(self, node s1, node t1, node s2, node t2):
		"""
		Changes the edge (s1, t1) into (s1, t2) and the edge (s2, t2) into (s2, t1).

		If there are edge weights or edge ids, they are preserved. Note that no check is performed if the swap is actually possible, i.e. does not generate duplicate edges.

		Parameters:
		-----------
		s1 : node
			Source node of the first edge
		t1 : node
			Target node of the first edge
		s2 : node
			Source node of the second edge
		t2 : node
			Target node of the second edge
		"""
		self._this.swapEdge(s1, t1, s2, t2)
		return self

	def compactEdges(self):
		"""
		Compact the edge storage, this should be called after executing many edge deletions.
		"""
		self._this.compactEdges()

	def sortEdges(self):
		"""
		Sorts the adjacency arrays by node id. While the running time is linear this
		temporarily duplicates the memory.
		"""
		self._this.sortEdges()

	def hasEdge(self, u, v):
		""" Checks if undirected edge {`u`,`v`} exists in the graph.

		Parameters:
		-----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.

		Returns:
		--------
		bool
			True if the edge exists, False otherwise.
		"""
		return self._this.hasEdge(u, v)

	def weight(self, u, v):
		""" Get edge weight of edge {`u` , `v`}. Returns 0 if edge does not exist.

		Parameters:
		-----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.

		Returns:
		--------
		edgeweight
			Edge weight of edge {`u` , `v`} or 0 if edge does not exist.
		"""
		return self._this.weight(u, v)

	def forNodes(self, object callback):
		""" Experimental node iterator interface

		Parameters:
		-----------
		callback : object
			Any callable object that takes the parameter node
		"""
		cdef NodeCallbackWrapper* wrapper
		try:
			wrapper = new NodeCallbackWrapper(callback)
			self._this.forNodes[NodeCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forNodesInRandomOrder(self, object callback):
		""" Experimental node iterator interface

		Parameters:
		-----------
		callback : object
			Any callable object that takes the parameter node
		"""
		cdef NodeCallbackWrapper* wrapper
		try:
			wrapper = new NodeCallbackWrapper(callback)
			self._this.forNodesInRandomOrder[NodeCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forNodePairs(self, object callback):
		""" Experimental node pair iterator interface

		Parameters:
		-----------
		callback : object
			Any callable object that takes the parameters (node, node)
		"""
		cdef NodePairCallbackWrapper* wrapper
		try:
			wrapper = new NodePairCallbackWrapper(callback)
			self._this.forNodePairs[NodePairCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forEdges(self, object callback):
		""" Experimental edge iterator interface

		Parameters:
		-----------
		callback : object
			Any callable object that takes the parameter (node, node, edgeweight, edgeid)
		"""
		cdef EdgeCallBackWrapper* wrapper
		try:
			wrapper = new EdgeCallBackWrapper(callback)
			self._this.forEdges[EdgeCallBackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forEdgesOf(self, node u, object callback):
		""" Experimental incident (outgoing) edge iterator interface

		Parameters:
		-----------
		u : node
			The node of which incident edges shall be passed to the callback
		callback : object
			Any callable object that takes the parameter (node, node, edgeweight, edgeid)
		"""
		cdef EdgeCallBackWrapper* wrapper
		try:
			wrapper = new EdgeCallBackWrapper(callback)
			self._this.forEdgesOf[EdgeCallBackWrapper](u, dereference(wrapper))
		finally:
			del wrapper

	def forInEdgesOf(self, node u, object callback):
		""" Experimental incident incoming edge iterator interface

		Parameters:
		-----------
		u : node
			The node of which incident edges shall be passed to the callback
		callback : object
			Any callable object that takes the parameter (node, node, edgeweight, edgeid)
		"""
		cdef EdgeCallBackWrapper* wrapper
		try:
			wrapper = new EdgeCallBackWrapper(callback)
			self._this.forInEdgesOf[EdgeCallBackWrapper](u, dereference(wrapper))
		finally:
			del wrapper

	def isWeighted(self):
		"""
		Returns:
		--------
		bool
			True if this graph supports edge weights other than 1.0.
		"""
		return self._this.isWeighted()

	def isDirected(self):
		return self._this.isDirected()

	def totalEdgeWeight(self):
		""" Get the sum of all edge weights.

		Returns:
		--------
		edgeweight
			The sum of all edge weights.
		"""
		return self._this.totalEdgeWeight()

	def numberOfSelfLoops(self):
		""" Get number of self-loops, i.e. edges {v, v}.
		Returns:
		--------
		count
			number of self-loops.
		"""
		return self._this.numberOfSelfLoops()

	def checkConsistency(self):
		"""
		Check for invalid graph states, such as multi-edges.
		"""
		return self._this.checkConsistency()

	def iterNodes(self):
		"""
		Iterates over the nodes of the graph.
		"""
		it = self._this.nodeRange().begin()
		while it != self._this.nodeRange().end():
			yield dereference(it)
			preincrement(it)

	def iterEdges(self):
		"""
		Iterates over the edges of the graph.
		"""
		it = self._this.edgeRange().begin()
		while it != self._this.edgeRange().end():
			yield dereference(it).u, dereference(it).v
			preincrement(it)

	def iterEdgesWeights(self):
		"""
		Iterates over the edges of the graph and their weights.
		"""
		it = self._this.edgeWeightRange().begin()
		while it != self._this.edgeWeightRange().end():
			yield dereference(it).u, dereference(it).v, dereference(it).weight
			preincrement(it)

	def iterNeighbors(self, u):
		"""
		Iterates over a range of the neighbors of a node.

		Parameters:
		-----------
		u : Node
		"""
		it = self._this.neighborRange(u).begin()
		while it != self._this.neighborRange(u).end():
			yield dereference(it)
			preincrement(it)

	def iterInNeighbors(self, u):
		"""
		Iterates over a range of the in-neighbors of a node.

		Parameters:
		-----------
		u : Node
		"""
		it = self._this.inNeighborRange(u).begin()
		while it != self._this.inNeighborRange(u).end():
			yield dereference(it)
			preincrement(it)

	def iterNeighborsWeights(self, u):
		"""
		Iterates over a range of the neighbors of a node including the edge weights.
		The iterator is not safe to use with unweighted graphs. To avoid unsafe behavior
		a runtime error will be thrown.

		Parameters:
		-----------
		u : Node
		"""
		if not self._this.isWeighted():
			raise RuntimeError("iterNeighborsWeights: Use this iterator only on weighted graphs.")
		
		it = self._this.weightNeighborRange(u).begin()
		while it != self._this.weightNeighborRange(u).end():
			yield dereference(it)
			preincrement(it)		

	def iterInNeighborsWeights(self, u):
		"""
		Iterates over a range of the in-neighbors of a node including the edge weights.
		The iterator is not safe to use with unweighted graphs. To avoid unsafe behavior
		a runtime error will be thrown.

		Parameters:
		-----------
		u : Node
		"""
		if not self._this.isWeighted():
			raise RuntimeError("iterInNeighborsWeights: Use this iterator only on weighted graphs.")

		it = self._this.weightInNeighborRange(u).begin()
		while it != self._this.weightInNeighborRange(u).end():
			yield dereference(it)
			preincrement(it)

cdef cppclass EdgeCallBackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, node v, edgeweight w, edgeid eid):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(u, v, w, eid)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef cppclass NodeCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(u)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef cppclass NodeDistCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, count dist):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(u, dist)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef cppclass NodePairCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, node v):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(u, v)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef class SpanningForest:
	""" Generates a spanning forest for a given graph

		Parameters:
		-----------
		G : networkit.Graph
			The graph.
		nodes : list
			A subset of nodes of `G` which induce the subgraph.
	"""
	cdef _SpanningForest* _this
	cdef Graph _G

	def __cinit__(self, Graph G not None):
		self._G = G
		self._this = new _SpanningForest(G._this)


	def __dealloc__(self):
		del self._this

	def run(self):
		"""
		Executes the algorithm.

		Returns:
		--------
		Algorithm:
			self
		"""
		self._this.run()
		return self

	def getForest(self):
		"""
		Returns the spanning forest.

		Returns:
		--------
		networkit.Graph
			The computed spanning forest
		"""
		return Graph().setThis(self._this.getForest())

cdef class RandomMaximumSpanningForest(Algorithm):
	"""
	Computes a random maximum-weight spanning forest using Kruskal's algorithm by randomizing the order of edges of the same weight.
	Parameters:
	-----------
	G : networkit.Graph
		The input graph.
	attribute : list
		If given, this edge attribute is used instead of the edge weights.
	"""

	def __cinit__(self, Graph G not None, vector[double] attribute = vector[double]()):
		self._G = G
		if attribute.empty():
			self._this = new _RandomMaximumSpanningForest(G._this)
		else:
			self._attribute = move(attribute)
			self._this = new _RandomMaximumSpanningForest(G._this, self._attribute)

	def getMSF(self, bool_t move):
		"""
		Gets the calculated maximum-weight spanning forest as graph.
		Parameters:
		-----------
		move : bool
			If the graph shall be moved out of the algorithm instance.

		Returns:
		--------
		networkit.Graph
			The calculated maximum-weight spanning forest.
		"""
		return Graph().setThis((<_RandomMaximumSpanningForest*>(self._this)).getMSF(move))

	def getAttribute(self, bool_t move = False):
		"""
		Get a bool attribute that indicates for each edge if it is part of the calculated maximum-weight spanning forest.
		This attribute is only calculated and can thus only be request if the supplied graph has edge ids.
		Parameters:
		-----------
		move : bool
			If the attribute shall be moved out of the algorithm instance.

		Returns:
		--------
		list
			The list with the bool attribute for each edge.
		"""
		return (<_RandomMaximumSpanningForest*>(self._this)).getAttribute(move)

	def inMSF(self, node u, node v = _none):
		"""
		Checks if the edge (u, v) or the edge with id u is part of the calculated maximum-weight spanning forest.
		Parameters:
		-----------
		u : node or edgeid
			The first node of the edge to check or the edge id of the edge to check
		v : node
			The second node of the edge to check (only if u is not an edge id)

		Returns:
		--------
		bool
			If the edge is part of the calculated maximum-weight spanning forest.
		"""
		if v == _none:
			return (<_RandomMaximumSpanningForest*>(self._this)).inMSF(u)
		else:
			return (<_RandomMaximumSpanningForest*>(self._this)).inMSF(u, v)

cdef class UnionMaximumSpanningForest(Algorithm):
	"""
	Union maximum-weight spanning forest algorithm, computes the union of all maximum-weight spanning forests using Kruskal's algorithm.

	Parameters:
	-----------
	G : networkit.Graph
		The input graph.
	attribute : list
		If given, this edge attribute is used instead of the edge weights.
	"""

	def __cinit__(self, Graph G not None, vector[double] attribute = vector[double]()):
		self._G = G

		if attribute.empty():
			self._this = new _UnionMaximumSpanningForest(G._this)
		else:
			self._this = new _UnionMaximumSpanningForest(G._this, attribute)

	def getUMSF(self, bool_t move = False):
		"""
		Gets the union of all maximum-weight spanning forests as graph.

		Parameters:
		-----------
		move : bool
			If the graph shall be moved out of the algorithm instance.

		Returns:
		--------
		networkit.Graph
			The calculated union of all maximum-weight spanning forests.
		"""
		return Graph().setThis((<_UnionMaximumSpanningForest*>(self._this)).getUMSF(move))

	def getAttribute(self, bool_t move = False):
		"""
		Get a bool attribute that indicates for each edge if it is part of any maximum-weight spanning forest.

		This attribute is only calculated and can thus only be request if the supplied graph has edge ids.

		Parameters:
		-----------
		move : bool
			If the attribute shall be moved out of the algorithm instance.

		Returns:
		--------
		list
			The list with the bool attribute for each edge.
		"""
		return (<_UnionMaximumSpanningForest*>(self._this)).getAttribute(move)

	def inUMST(self, node u, node v = _none):
		"""
		Checks if the edge (u, v) or the edge with id u is part of any maximum-weight spanning forest.

		Parameters:
		-----------
		u : node or edgeid
			The first node of the edge to check or the edge id of the edge to check
		v : node
			The second node of the edge to check (only if u is not an edge id)

		Returns:
		--------
		bool
			If the edge is part of any maximum-weight spanning forest.
		"""
		if v == _none:
			return (<_UnionMaximumSpanningForest*>(self._this)).inUMSF(u)
		else:
			return (<_UnionMaximumSpanningForest*>(self._this)).inUMSF(u, v)
