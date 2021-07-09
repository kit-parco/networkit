# distutils: language=c++

from libc.stdint cimport uint64_t
from libc.stdint cimport uint8_t

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.map cimport map
from libcpp.unordered_map cimport unordered_map

import os
import logging
import numpy
import scipy.io
import fnmatch
from enum import Enum

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node

from .dynamics import DGSWriter, DGSStreamParser
from .graph cimport _Graph, Graph
from .graph import Graph as __Graph
from .structures cimport _Cover, Cover, _Partition, Partition
from .GraphMLIO import GraphMLReader, GraphMLWriter
from .GEXFIO import GEXFReader, GEXFWriter
from . import algebraic
from .support import MissingDependencyError

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Graph move( _Graph t ) nogil # specialized declaration as general declaration disables template argument deduction and doesn't work
	_Partition move( _Partition t) nogil

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes

cdef extern from "<networkit/io/GraphReader.hpp>":

	cdef cppclass _GraphReader "NetworKit::GraphReader":
		_GraphReader() nogil except +
		_Graph read(string path) nogil except +

cdef class GraphReader:
	""" Abstract base class for graph readers"""

	cdef _GraphReader* _this

	def __init__(self, *args, **kwargs):
		if type(self) == GraphReader:
			raise RuntimeError("Error, you may not use GraphReader directly, use a sub-class instead")

	def __cinit__(self, *args, **kwargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

	def read(self, path):
		cdef string cpath = stdstring(path)
		cdef _Graph result

		with nogil:
			result = move(self._this.read(cpath)) # extra move in order to avoid copying the internal variable that is used by Cython
		return Graph(0).setThis(result)

cdef extern from "<networkit/io/GraphReader.hpp>" namespace "NetworKit::GraphReader":

	cdef enum _MultipleEdgesHandling "NetworKit::GraphReader::MultipleEdgesHandling":
		DISCARD_EDGES,
		SUM_WEIGHTS_UP,
		KEEP_MINIMUM_WEIGHT

class MultipleEdgesHandling:
	DiscardEdges = DISCARD_EDGES
	SumWeightsUp = SUM_WEIGHTS_UP
	KeepMinimumWeight = KEEP_MINIMUM_WEIGHT

cdef extern from "<networkit/io/GraphWriter.hpp>":

	cdef cppclass _GraphWriter "NetworKit::GraphWriter":
		_GraphWriter() nogil except +
		void write(_Graph G, string path) nogil except +


cdef class GraphWriter:
	"""
	Abstract base class for graph writers
	"""

	cdef _GraphWriter *_this

	def __init__(self, *args, **kwargs):
		if type(self) == GraphWriter:
			raise RuntimeError("Error, you may not use GraphWriter directly, use a sub-class instead")

	def __cinit__(self, *args, **kwargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

	def write(self, Graph G not None, path):
		"""
		Write the graph to a file.

		Parameters:
		-----------
		G     : networkit.Graph
			The graph to write
		paths : str
			The output path
		"""
		assert path != None
		cdef string c_path = stdstring(path)
		with nogil:
			self._this.write(G._this, c_path)
		return self

cdef extern from "<networkit/io/METISGraphReader.hpp>":

	cdef cppclass _METISGraphReader "NetworKit::METISGraphReader" (_GraphReader):
		_METISGraphReader() nogil except +

cdef class METISGraphReader(GraphReader):
	""" Reads the METIS adjacency file format [1]. If the Fast reader fails,
		use readGraph(path, graphio.formats.metis) as an alternative.
		[1]: http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
	"""
	def __cinit__(self):
		self._this = new _METISGraphReader()

cdef extern from "<networkit/io/NetworkitBinaryReader.hpp>":
	cdef cppclass _NetworkitBinaryReader "NetworKit::NetworkitBinaryReader" (_GraphReader):
		_NetworkitBinaryReader() except +
		_Graph readState(string state) except +

cdef class NetworkitBinaryReader(GraphReader):
	"""
	Reads a graph written in the custom Networkit format documented in cpp/io/NetworkitGraph.md
	"""

	def __cinit__(self):
		self._this = new _NetworkitBinaryReader()
	
	def readState(self, state):
		return Graph().setThis((<_NetworkitBinaryReader*>(self._this)).readState(state))

cdef extern from "<networkit/io/NetworkitBinaryWriter.hpp>":
	cdef cppclass _NetworkitBinaryWriter "NetworKit::NetworkitBinaryWriter" (_GraphWriter):
		_NetworkitBinaryWriter() except +
		string writeState(_Graph G) except +

cdef class NetworkitBinaryWriter(GraphWriter):
	def __cinit__(self):
		self._this = new _NetworkitBinaryWriter()
	
	def writeState(self, Graph G not None):
		return (<_NetworkitBinaryWriter*>(self._this)).writeState(G._this)

cdef extern from "<networkit/io/GraphToolBinaryReader.hpp>":

	cdef cppclass _GraphToolBinaryReader "NetworKit::GraphToolBinaryReader" (_GraphReader):
		_GraphToolBinaryReader() except +

cdef class GraphToolBinaryReader(GraphReader):
	""" Reads the binary file format defined by graph-tool[1].
		[1]: http://graph-tool.skewed.de/static/doc/gt_format.html
	"""
	def __cinit__(self):
		self._this = new _GraphToolBinaryReader()

cdef extern from "<networkit/io/ThrillGraphBinaryReader.hpp>":

	cdef cppclass _ThrillGraphBinaryReader "NetworKit::ThrillGraphBinaryReader" (_GraphReader):
		_ThrillGraphBinaryReader(count n) except +
		_Graph read(vector[string] paths) nogil except +

cdef class ThrillGraphBinaryReader(GraphReader):
	"""
	Reads a graph format consisting of a serialized DIA of vector<uint32_t> from thrill.
	When the number of nodes is given, reading the graph is more efficient.
	Otherwise nodes are added to the graph as they are encountered.
	Edges must be present only in one direction.

	Parameters:
	-----------
	n : count
		The number of nodes
	"""
	def __cinit__(self, count n = 0):
		self._this = new _ThrillGraphBinaryReader(n)

	"""
	Read the graph from one or multiple files

	Parameters:
	-----------
	paths : str or list[str]
		The input path(s)
	"""
	def read(self, paths):
		cdef vector[string] c_paths

		if isinstance(paths, str):
			c_paths.push_back(stdstring(paths))
		else:
			c_paths.reserve(len(paths))

			for p in paths:
				c_paths.push_back(stdstring(p))

		cdef _Graph result

		with nogil:
			result = move((<_ThrillGraphBinaryReader*>(self._this)).read(c_paths)) # extra move in order to avoid copying the internal variable that is used by Cython

		return Graph(0).setThis(result)

cdef extern from "<networkit/io/ThrillGraphBinaryWriter.hpp>":

	cdef cppclass _ThrillGraphBinaryWriter "NetworKit::ThrillGraphBinaryWriter" (_GraphWriter):
		_ThrillGraphBinaryWriter() except +

cdef class ThrillGraphBinaryWriter(GraphWriter):
	"""
	Writes a graph format consisting of a serialized DIA of vector<uint32_t> from Thrill.
	Edges are written only in one direction.
	"""

	def __cinit__(self):
		self._this = new _ThrillGraphBinaryWriter()

cdef extern from "<networkit/io/EdgeListReader.hpp>":

	cdef cppclass _EdgeListReader "NetworKit::EdgeListReader"(_GraphReader):
		_EdgeListReader() except +
		_EdgeListReader(char separator, node firstNode, string commentPrefix, bool_t continuous, bool_t directed)
		map[string,node] getNodeMap() except +

cdef class EdgeListReader(GraphReader):
	""" Reads a graph from various text-based edge list formats.

	EdgeListReader(self, separator, firstNode, commentPrefix="#", continuous=True, directed=False)

	A line has to contain two or three entries separated with the separator symbol (one ASCII character).
	If at least one line contains three entries, the generated graph will be weighted and
	each line with only two fields will be interpreted as weight 1.0.

	A file may contain the same edge multiple times; then, the weight of the first
	occurrence is used.

	Undirected graphs need to include an edge only in one direction, i.e. edge {u, v} may
	be represented by (u, v) or (v, u) or both (again, only the first occurrence is used).

	If the input file contains non-continuous node ids with large gaps or non-integer node labels,
	set the parameter continuous to False. Then, gaps are automatically removed and node ids are
	reassigned to [0, n) where n is the number of nodes in the graph. The mapping will be arbitrary
	and can be accessed using getNodeMap().

	To shift continuous integer node labels which are not zero-indexed, set firstNode to
	the smallest id used in the file. firstNode will be ignored in the non-continuous case.

	The file may also include line comments which start with the commentPrefix.

	Parameters:
	-----------
	separator : char
		The separator character. Must have length of exactly one.
	firstNode : node
		The id of the first node, this value will be subtracted from all node ids
	commentPrefix : string
		Lines starting with this prefix will be ignored
	continuous : bool
		File uses continuous node ids.
	directed : bool
		Treat input file as a directed graph.
	"""
	def __cinit__(self, separator, firstNode, commentPrefix="#", continuous=True, directed=False):
		if len(separator) != 1 or ord(separator[0]) > 255:
			raise RuntimeError("separator has to be exactly one ascii character");

		self._this = new _EdgeListReader(stdstring(separator)[0], firstNode, stdstring(commentPrefix), continuous, directed)

	def getNodeMap(self):
		""" Returns mapping of non-continuous files.

		The mapping is returned as dict (string -> node) projecting the original
		labels (as strings) to the reassigned integer node ids.
		"""
		cdef map[string,node] cResult = (<_EdgeListReader*>(self._this)).getNodeMap()
		result = dict()
		for elem in cResult:
			result[(elem.first).decode("utf-8")] = elem.second
		return result

cdef extern from "<networkit/io/KONECTGraphReader.hpp>":

	cdef cppclass _KONECTGraphReader "NetworKit::KONECTGraphReader"(_GraphReader):
		_KONECTGraphReader() except +
		_KONECTGraphReader(bool_t remapNodes, _MultipleEdgesHandling handlingmethod)

cdef class KONECTGraphReader(GraphReader):
	""" Reader for the KONECT graph format, which is described in detail on the KONECT website[1].

		[1]: http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
	"""
	def __cinit__(self, remapNodes = False, handlingmethod = MultipleEdgesHandling.DiscardEdges):
		self._this = new _KONECTGraphReader(remapNodes, handlingmethod)

cdef extern from "<networkit/io/GMLGraphReader.hpp>":

	cdef cppclass _GMLGraphReader "NetworKit::GMLGraphReader"(_GraphReader):
		_GMLGraphReader() except +

cdef class GMLGraphReader(GraphReader):
	""" Reader for the GML graph format, which is documented here [1].

		[1]: http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf
 	"""
	def __cinit__(self):
		self._this = new _GMLGraphReader()

cdef extern from "<networkit/io/METISGraphWriter.hpp>":

	cdef cppclass _METISGraphWriter "NetworKit::METISGraphWriter" (_GraphWriter):
		_METISGraphWriter() except +


cdef class METISGraphWriter(GraphWriter):
	""" Writes graphs in the METIS format"""

	def __cinit__(self):
		self._this = new _METISGraphWriter()
cdef extern from "<networkit/io/GraphToolBinaryWriter.hpp>":

	cdef cppclass _GraphToolBinaryWriter "NetworKit::GraphToolBinaryWriter" (_GraphWriter):
		_GraphToolBinaryWriter() except +

cdef class GraphToolBinaryWriter(GraphWriter):
	""" Reads the binary file format defined by graph-tool[1].
		[1]: http://graph-tool.skewed.de/static/doc/gt_format.html
	"""
	def __cinit__(self):
		self._this = new _GraphToolBinaryWriter()

cdef extern from "<networkit/io/DotGraphWriter.hpp>":

	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter" (_GraphWriter):
		_DotGraphWriter() except +

cdef class DotGraphWriter(GraphWriter):
	""" Writes graphs in the .dot/GraphViz format"""
	def __cinit__(self):
		self._this = new _DotGraphWriter()

cdef extern from "<networkit/io/GMLGraphWriter.hpp>":

	cdef cppclass _GMLGraphWriter "NetworKit::GMLGraphWriter" (_GraphWriter):
		_GMLGraphWriter() except +


cdef class GMLGraphWriter(GraphWriter):
	""" Writes a graph and its coordinates as a GML file.[1]
		[1] http://svn.bigcat.unimaas.nl/pvplugins/GML/trunk/docs/gml-technical-report.pdf """

	def __cinit__(self):
		self._this = new _GMLGraphWriter()

cdef extern from "<networkit/io/EdgeListWriter.hpp>":

	cdef cppclass _EdgeListWriter "NetworKit::EdgeListWriter" (_GraphWriter):
		_EdgeListWriter() except +
		_EdgeListWriter(char separator, node firstNode, bool_t bothDirections) except +

cdef class EdgeListWriter(GraphWriter):
	""" Writes graphs in various edge list formats.

	Parameters:
	-----------
	separator : string
		The separator character.
	firstNode : node
		The id of the first node, this value will be added to all node ids
	bothDirections : bool, optional
		If undirected edges shall be written in both directions, i.e., as symmetric directed graph (default: false)
	"""

	def __cinit__(self, separator, firstNode, bool_t bothDirections = False):
		cdef char sep = stdstring(separator)[0]
		self._this = new _EdgeListWriter(sep, firstNode, bothDirections)

cdef extern from "<networkit/io/LineFileReader.hpp>":

	cdef cppclass _LineFileReader "NetworKit::LineFileReader":
		_LineFileReader() except +
		vector[string] read(string path)


cdef class LineFileReader:
	""" Reads a file and puts each line in a list of strings """
	cdef _LineFileReader _this

	def read(self, path):
		return self._this.read(stdstring(path))


cdef extern from "<networkit/io/SNAPGraphWriter.hpp>":
	cdef cppclass _SNAPGraphWriter "NetworKit::SNAPGraphWriter" (_GraphWriter):
		_SNAPGraphWriter() except +

cdef class SNAPGraphWriter(GraphWriter):
	""" Writes graphs in a format suitable for the Georgia Tech SNAP software [1]
		[1]: http://snap-graph.sourceforge.net/
	"""

	def __cinit__(self):
		self._this = new _SNAPGraphWriter()

cdef extern from "<networkit/io/SNAPGraphReader.hpp>":

	cdef cppclass _SNAPGraphReader "NetworKit::SNAPGraphReader"(_GraphReader):
		_SNAPGraphReader() except +
		_SNAPGraphReader(bool_t directed, bool_t remapNodes, count nodeCount)

cdef class SNAPGraphReader(GraphReader):
	""" Reads a graph from the SNAP graph data collection [1]
		[1]: http://snap.stanford.edu/data/index.html
	"""
	def __cinit__(self, directed = False, remapNodes = True, nodeCount = 0):
		self._this = new _SNAPGraphReader(directed, remapNodes, nodeCount)



cdef extern from "<networkit/io/PartitionReader.hpp>":

	cdef cppclass _PartitionReader "NetworKit::PartitionReader":
		_PartitionReader() except +
		_Partition read(string path) except +


cdef class PartitionReader:
	""" Reads a partition from a file.
		File format: line i contains subset id of element i.
	 """
	cdef _PartitionReader _this

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))


cdef extern from "<networkit/io/PartitionWriter.hpp>":

	cdef cppclass _PartitionWriter "NetworKit::PartitionWriter":
		_PartitionWriter() except +
		void write(_Partition, string path) nogil except +


cdef class PartitionWriter:
	""" Writes a partition to a file.
		File format: line i contains subset id of element i.
	 """
	cdef _PartitionWriter _this

	def write(self, Partition zeta, path):
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(zeta._this, cpath)

cdef extern from "<networkit/io/BinaryPartitionReader.hpp>":

	cdef cppclass _BinaryPartitionReader "NetworKit::BinaryPartitionReader":
		_BinaryPartitionReader() except +
		_BinaryPartitionReader(uint8_t width) except +
		_Partition read(string path) except +


cdef class BinaryPartitionReader:
	"""
	Reads a partition from a binary file that contains an unsigned integer
	of the given width for each node.

	Parameters:
	-----------
	width : int
		the width of the unsigned integer in bytes (4 or 8)

	"""
	cdef _BinaryPartitionReader _this

	def __cinit__(self, uint8_t width=4):
		self._this = _BinaryPartitionReader(width)

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))

cdef extern from "<networkit/io/BinaryPartitionWriter.hpp>":

	cdef cppclass _BinaryPartitionWriter "NetworKit::BinaryPartitionWriter":
		_BinaryPartitionWriter() except +
		_BinaryPartitionWriter(uint8_t width) except +
		_Partition write(_Partition zeta, string path) nogil except +

cdef class BinaryPartitionWriter:
	"""
	Writes a partition to a file to contains a binary list of partition ids.
	Partition ids are unsigned integers.

        Parameters:
	-----------
	width : int
		the width of the unsigned integer in bytes (4 or 8)

	"""
	cdef _BinaryPartitionWriter _this

	def __cinit__(self, uint8_t width=4):
		self._this = _BinaryPartitionWriter(width)

	def write(self, Partition P not None, path):
		"""
		Write the partition to the given file.

		Parameters:
		-----------
		path : str
			The output path
		"""
		cdef string c_path = stdstring(path)

		with nogil:
			self._this.write(P._this, c_path)

		return self

cdef extern from "<networkit/io/EdgeListPartitionReader.hpp>":

	cdef cppclass _EdgeListPartitionReader "NetworKit::EdgeListPartitionReader":
		_EdgeListPartitionReader() except +
		_EdgeListPartitionReader(node firstNode, char sepChar) except +
		_Partition read(string path) except +


cdef class EdgeListPartitionReader:
	""" Reads a partition from an edge list type of file
	 """
	cdef _EdgeListPartitionReader _this

	def __cinit__(self, node firstNode=1, sepChar = '\t'):
		self._this = _EdgeListPartitionReader(firstNode, stdstring(sepChar)[0])

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))

cdef extern from "<networkit/io/BinaryEdgeListPartitionReader.hpp>":

	cdef cppclass _BinaryEdgeListPartitionReader "NetworKit::BinaryEdgeListPartitionReader":
		_BinaryEdgeListPartitionReader() except +
		_BinaryEdgeListPartitionReader(node firstNode, uint8_t width) except +
		_Partition read(string path) nogil except +
		_Partition read(vector[string] paths) nogil except +

cdef class BinaryEdgeListPartitionReader:
	"""
	Reads a partition file that contains a binary list of pairs (node, partition(node)).
	It is assumed that all integers are unsigned.

	Parameters:
	-----------
	firstNode : node
		The id of the first node, this is subtracted from all read node ids
	width : int
		The width of the unsigned integer in bytes (4 or 8)
	"""
	cdef _BinaryEdgeListPartitionReader _this

	def __cinit__(self, node firstNode=0, uint8_t width=4):
		self._this = _BinaryEdgeListPartitionReader(firstNode, width)

	def read(self, paths):
		"""
		Read the partition from one or multiple files

		Parameters:
		-----------
		paths : str or list[str]
			The input path(s)
		"""
		cdef vector[string] c_paths

		if isinstance(paths, str):
			c_paths.push_back(stdstring(paths))
		else:
			c_paths.reserve(len(paths))

			for p in paths:
				c_paths.push_back(stdstring(p))

		cdef _Partition result

		with nogil:
			result = move(self._this.read(c_paths)) # extra move in order to avoid copying the internal variable that is used by Cython

		return Partition().setThis(result)

cdef extern from "<networkit/io/BinaryEdgeListPartitionWriter.hpp>":

	cdef cppclass _BinaryEdgeListPartitionWriter "NetworKit::BinaryEdgeListPartitionWriter":
		_BinaryEdgeListPartitionWriter() except +
		_BinaryEdgeListPartitionWriter(node firstNode, uint8_t width) except +
		_Partition write(_Partition P, string path) nogil except +

cdef class BinaryEdgeListPartitionWriter:
	"""
	Writes a partition file that contains a binary list of pairs (node, partition(node)).

	Parameters:
	-----------
	firstNode : node
		The id of the first node, this is added to all writen node ids
	width : int
		The width of the unsigned integer in bytes (4 or 8)
	"""
	cdef _BinaryEdgeListPartitionWriter _this

	def __cinit__(self, node firstNode=0, uint8_t width=4):
		self._this = _BinaryEdgeListPartitionWriter(firstNode, width)

	def write(self, Partition P not None, path):
		"""
		Write the partition to the given file.

		Parameters:
		-----------
		path : str
			The output path
		"""
		cdef string c_path = stdstring(path)

		with nogil:
			self._this.write(P._this, c_path)

		return self

cdef extern from "<networkit/io/SNAPEdgeListPartitionReader.hpp>":

	cdef cppclass _SNAPEdgeListPartitionReader "NetworKit::SNAPEdgeListPartitionReader":
		_SNAPEdgeListPartitionReader() except +
		_Cover read(string path, unordered_map[node,node] nodeMap,_Graph G) except +
#		_Partition readWithInfo(string path, count nNodes) except +

cdef class SNAPEdgeListPartitionReader:
	""" Reads a partition from a SNAP 'community with ground truth' file
	 """
	cdef _SNAPEdgeListPartitionReader _this

	def read(self,path, nodeMap, Graph G):
		cdef unordered_map[node,node] cNodeMap
		for (key,val) in nodeMap:
			cNodeMap[key] = val
		return Cover().setThis(self._this.read(stdstring(path), cNodeMap, G._this))

cdef extern from "<networkit/io/CoverReader.hpp>":

	cdef cppclass _CoverReader "NetworKit::CoverReader":
		_CoverReader() except +
		_Cover read(string path,_Graph G) except +

cdef class CoverReader:
	""" Reads a cover from a file
		File format: each line contains the space-separated node ids of a community
	 """
	cdef _CoverReader _this

	def read(self, path, Graph G):
		return Cover().setThis(self._this.read(stdstring(path), G._this))

cdef extern from "<networkit/io/CoverWriter.hpp>":

	cdef cppclass _CoverWriter "NetworKit::CoverWriter":
		_CoverWriter() except +
		void write(_Cover, string path) nogil except +


cdef class CoverWriter:
	""" Writes a partition to a file.
		File format: each line contains the space-separated node ids of a community
	 """
	cdef _CoverWriter _this

	def write(self, Cover zeta, path):
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(zeta._this, cpath)

cdef extern from "<networkit/io/EdgeListCoverReader.hpp>":

	cdef cppclass _EdgeListCoverReader "NetworKit::EdgeListCoverReader":
		_EdgeListCoverReader() except +
		_EdgeListCoverReader(node firstNode) except +
		_Cover read(string path, _Graph G) except +


cdef class EdgeListCoverReader:
	""" Reads a cover from an edge list type of file
		File format: each line starts with a node id and continues with a list of the communities the node belongs to
	 """
	cdef _EdgeListCoverReader _this

	def __cinit__(self, firstNode=1):
		self._this = _EdgeListCoverReader(firstNode)

	def read(self, path, Graph G):
		return Cover().setThis(self._this.read(stdstring(path), G._this))

class __AutoNumber(Enum):
	def __new__(cls):
		value = len(cls.__members__) + 1
		obj = object.__new__(cls)
		obj._value_ = value
		return obj


class Format(__AutoNumber):
	""" Simple enumeration class to list supported file types. Currently supported
	file types: SNAP, EdgeListSpaceZero, EdgeListSpaceOne, EdgeListTabZero, EdgeListTabOne,
	METIS, GraphML, GEXF, GML, EdgeListCommaOne, GraphViz, DOT, EdgeList, LFR, KONEC, GraphToolBinary,
			NetworkitBinary"""
	SNAP = ()
	EdgeListSpaceZero = ()
	EdgeListSpaceOne = ()
	EdgeListTabZero = ()
	EdgeListTabOne = ()
	METIS = ()
	GraphML = ()
	GEXF = ()
	GML = ()
	EdgeListCommaOne = ()
	GraphViz = ()
	DOT = ()
	EdgeList = ()
	LFR = ()
	KONECT = ()
	GraphToolBinary = ()
	MAT = ()
	ThrillBinary = ()
	NetworkitBinary = ()

# reading

def getReader(fileformat, *kargs, **kwargs):
	#define your [edgelist] reader here:
	readers = {
		Format.METIS:			METISGraphReader(),
		Format.GraphML:			GraphMLReader(),
		Format.GEXF:			GEXFReader(),
		Format.SNAP:			SNAPGraphReader(),
		Format.EdgeListCommaOne:	EdgeListReader(',',1,),
		Format.EdgeListSpaceOne:	EdgeListReader(' ',1),
		Format.EdgeListSpaceZero:	EdgeListReader(' ',0),
		Format.EdgeListTabOne:		EdgeListReader('\t',1),
		Format.EdgeListTabZero:		EdgeListReader('\t',0),
		Format.LFR:			EdgeListReader('\t',1),
		Format.KONECT:			KONECTGraphReader(),
		Format.GML:			GMLGraphReader(),
		Format.GraphToolBinary:		GraphToolBinaryReader(),
		Format.MAT:			MatReader(),
		Format.ThrillBinary:		ThrillGraphBinaryReader(),
		Format.NetworkitBinary:         NetworkitBinaryReader()
	}

	# special case for custom Edge Lists
	if fileformat == Format.EdgeList:
		if "continuous" in kwargs and kwargs["continuous"] == False:
			kwargs["firstNode"] = 0
		reader = EdgeListReader(*kargs, **kwargs)

	else:
		if fileformat not in readers:
			raise Exception("unrecognized format/format not supported as input: {0}".format(fileformat))
		reader = readers[fileformat]#(**kwargs)

	return reader


def readGraph(path, fileformat, *kargs, **kwargs):
	""" Read graph file in various formats and return a NetworKit::Graph
	    Parameters:
		- fileformat: An element of the Format enumeration. Currently supported file types:
		SNAP, EdgeListSpaceZero, EdgeListSpaceOne, EdgeListTabZero, EdgeListTabOne, METIS,
		GraphML, GEXF, GML, EdgeListCommaOne, GraphViz, DOT, EdgeList, LFR, KONECT, GraphToolBinary, ThrillBinary
		- **kwargs: in case of a custom edge list, pass the genereic Fromat.EdgeList accompanied by
			the defining paramaters as follows:
			"separator=CHAR, firstNode=NODE, commentPrefix=STRING, continuous=BOOL, directed=BOOL"
			commentPrefix='#', continuous=True and directed=False are optional because of their default values;
			firstNode is not needed when continuous=True.
	"""
	reader = getReader(fileformat, *kargs, **kwargs)

	if ("~" in path):
		path = os.path.expanduser(path)
		print("path expanded to: {0}".format(path))
	if not os.path.isfile(path):
		raise IOError("{0} is not a file".format(path))
	else:
		with open(path, "r") as file:    # catch a wrong path before it crashes the interpreter
			try:
				G = reader.read(path)
				return G
			except Exception as e:
				raise IOError("{0} is not a valid {1} file: {2}".format(path,fileformat,e))
	return None

def readGraphs(dirPath, pattern, fileformat, some=None, exclude=None, **kwargs):
	"""
	Read all graph files contained in a directory whose filename contains the pattern, return a dictionary of name to Graph object.
    Parameters:
	- pattern: unix-style string pattern
	- fileformat: An element of the Format enumeration
	- some: restrict number of graphs to be read
	- **kwargs: in case of a custom edge list, provide the defining paramaters as follows:
		"separator=CHAR, firstNode=NODE, commentPrefix=STRING, continuous=BOOL"
		commentPrefix and continuous are optional
	"""
	graphs = {}
	graph_id = 0
	for root, dirs, files in os.walk(dirPath):
		for file in files:
			if fnmatch.fnmatch(file, pattern):
				if (exclude is None) or (not fnmatch.fnmatch(file, exclude)):
					G = readGraph(os.path.join(root, file), fileformat, **kwargs)
					graphs[graph_id] = G
					graph_id += 1
					if some:
						if len(graphs) == some:
							return graphs
	return graphs


class MatReader:
	def __init__(self, key = "G"):
		self.key = key

	def read(self, path):
		return readMat(path, self.key)

def readMat(path, key="G"):
	""" Reads a Graph from a matlab object file containing an adjacency matrix and returns a NetworKit::Graph
		Parameters:
		- key: The key of the adjacency matrix in the matlab object file (default: A)"""
	matlabObject = scipy.io.loadmat(path)
	# result is a dictionary of variable names and objects, representing the matlab object
	if key in matlabObject:
		A = matlabObject[key]
	else:
		raise Exception("Key {0} not found in the matlab object file".format(key))
	(n, n2) = A.shape
	if n != n2:
		raise Exception("this ({0}x{1}) matrix is not square".format(n, n2))
#	if not numpy.array_equal(A, A.transpose): # FIXME this is slow and doesn't work as expected, seems to be False for valid inputs
#		logging.warning("the adjacency matrix is not symmetric")
	G = __Graph(n)
	nz = A.nonzero()
	for (u,v) in zip(nz[0], nz[1]):
		if not G.hasEdge(u, v):
			G.addEdge(u, v)
	return G

class MatWriter:
	def __init__(self, key="G"):
		self.key = key

	def write(self, G, path, key="G"):
		writeMat(G, path, key)

def writeMat(G, path, key="G"):
	""" Writes a NetworKit::Graph to a Matlab object file.
		Parameters:
		-----------
		- G: The graph
		- path: Path of the matlab file
		- key: Dictionary Key
	"""
	matrix = algebraic.adjacencyMatrix(G, matrixType='sparse')
	scipy.io.savemat(path, {key : matrix})


# writing
def getWriter(fileformat, *kargs, **kwargs):
	writers =	{
		Format.METIS:			METISGraphWriter(),
		Format.GraphML:			GraphMLWriter(),
		Format.GEXF:			GEXFWriter(),
		Format.SNAP:			SNAPGraphWriter(),
		Format.EdgeListCommaOne:	EdgeListWriter(',',1,),
		Format.EdgeListSpaceOne:	EdgeListWriter(' ',1),
		Format.EdgeListSpaceZero:	EdgeListWriter(' ',0),
		Format.EdgeListTabOne:		EdgeListWriter('\t',1),
		Format.EdgeListTabZero:		EdgeListWriter('\t',0),
		Format.GraphViz:		DotGraphWriter(),
		Format.DOT:			DotGraphWriter(),
		Format.GML:			GMLGraphWriter(),
		Format.LFR:			EdgeListWriter('\t',1),
		Format.GraphToolBinary:         GraphToolBinaryWriter(),
		Format.MAT:			MatWriter(),
		Format.ThrillBinary:		ThrillGraphBinaryWriter(),
		Format.NetworkitBinary:         NetworkitBinaryWriter()
	}

	# special case for custom Edge Lists
	if fileformat == Format.EdgeList:
		return EdgeListWriter(*kargs, **kwargs)

	else:
		if fileformat not in writers:
			raise Exception("format {0} currently not supported".format(fileformat))

		return writers[fileformat]#(**kwargs)

def writeGraph(G, path, fileformat, *kargs, **kwargs):
	""" Write graph to various output formats.

	Paramaters:
	-----------
	- G:			a graph
	- path: 		output path
	- fileformat: 	an element of the Format enumeration

	"""

	dirname = os.path.dirname(os.path.realpath(path))
	# the given file path does not exist yet
	if not os.path.isfile(path):
		# check write permissions on the directory
		if not os.access(dirname, os.W_OK):
			# we may not write on this directory, raise Error
			raise IOError("No permission to write")
		# else everthing is alright
	else:
		# the given path points to a file
		if not os.access(path, os.W_OK):
			raise IOError("No permission to write")
		else:
			logging.warning("overriding given file")
	writer = getWriter(fileformat, *kargs, **kwargs)
	writer.write(G, path)
	logging.info("wrote graph {0} to file {1}".format(G, path))


class GraphConverter:

	def __init__(self, reader, writer):
		self.reader = reader
		self.writer = writer

	def convert(self, inPath, outPath):
		G = self.reader.read(inPath)
		self.writer.write(G, outPath)

	def __str__(self):
		return "GraphConverter: {0} => {0}".format(self.reader, self.writer)

def getConverter(fromFormat, toFormat):
	reader = getReader(fromFormat)
	writer = getWriter(toFormat)
	return GraphConverter(reader, writer)


def convertGraph(fromFormat, toFormat, fromPath, toPath=None):
	converter = getConverter(fromFormat, toFormat)
	if toPath is None:
		toPath = "{0}.{1}.graph".format(fromPath.split(".")[0], toFormat)
	converter.convert(fromPath, toPath)
	print("converted {0} to {1}".format(fromPath, toPath))

# dynamic

def readStream(path, mapped=True, baseIndex=0):
	"""
		Read a graph event stream from a file.
	"""
	return DGSStreamParser(path, mapped, baseIndex).getStream()

def writeStream(stream, path):
	"""
		Write a graph event stream to a file.
	"""
	DGSWriter().write(stream, path)
