from libc.stdint cimport uint64_t
from libcpp.map cimport map
from libcpp.set cimport set

ctypedef uint64_t index
ctypedef index node

from .graph cimport _Graph, Graph

cdef extern from "<networkit/scd/SelectiveCommunityDetector.hpp>":

	cdef cppclass _SelectiveCommunityDetector "NetworKit::SelectiveCommunityDetector":
		_SelectiveCommunityDetector(_Graph G) except +
		map[node, set[node]] run(set[node] seeds) except +
		set[node] expandOneCommunity(node seed) except +
		set[node] expandOneCommunity(set[node] seeds) except +

cdef class SelectiveCommunityDetector:
	cdef _SelectiveCommunityDetector *_this
	cdef Graph _G
