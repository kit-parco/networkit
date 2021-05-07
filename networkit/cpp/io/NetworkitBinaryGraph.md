NetworkitBinaryGraph format
==========================
This is a documentation of the NetworkitBinaryGraph file format. 
The file starts with a header followed by the base data, adjacency lists and weights respectively. 
The rest of this documentation describes the aforementioned blocks in detail.

Header
---------
The file begins with a header that is structured as follows:
```
struct Header {

    char magic[8];
    uint64_t checksum;
    uint64_t features; 
    uint64_t nodes;
    uint64_t chunks;
    uint64_t offsetBaseData;
    uint64_t offsetAdjLists;
    uint64_t offsetTranspose;
    uint64_t offsetWeightLists;
    uint64_t offsetWeightTranspose;
};
```
- magic: A constant value used to identify the file format version.
    - The current version is '*nkbg002*' which supports weighted, undirected and directed graphs.
- checksum: Currently not used
- features: Contains the graph information bitwise
    - Bit 0 : directed or undirected
    - Bit 1-3 : the weight format of the graph. Below are the possible formats:
		- 0 = Graph is unweighted
        - 1 = Weights are unsigned integers
        - 2 = Weights are signed integers
        - 3 = Weights are doubles
        - 4 = Weights are floats
	- Bit 4 = proper edge indexing
- nodes : The number of nodes the graph has
- chunks: The number of chunks the nodes have been divided in
- offsetBaseData: Offset of base data in the file 
- offsetAdjLists: Offset of the adjacency lists in the file
- offsetTranspose: Offset of the transposed adjaceny lists in the file
- offsetWeightLists: Offset of the adjacency weights in the file
- offsetWeightTranspose: Offset of the transposed adjacency weights in the file

All offsets are relative to the beginning of the section.

Base data
------------
```
uint64_t nodeFlags[nodes]: Flags storing information about a node
uint64_t firstVertex[chunks-1]: The index of the first vertex of each chunk excluding the first chunk
uint64_t edgeIndex[edges]:
```
Adjacency lists
-----------------
```
uint64_t nrOfEdges: the total number of edges in the block
uint64_t offset[chunks-1]: Offset of the file where the adjacency list of 
the firstVertex of each chunk relative to data starts
varint data [...]: Varint encoded adjaceny lists  
```
Transpose lists
-----------------
```
uint64_t nrOfEdges: the total number of edges in the block
uint64_t offset[chunks-1]: Offset of the file where the tranpose list of 
the firstVertex of each chunk relative to data starts
varint data [...]: Varint encoded transpose lists  

```
Weight lists
--------------------
```
uint64_t offset[chunks-1]: Offset of the file where the weights are
Depending on the type of weights:
 - unsigned/signed weights: varint data [...]: Varint encoded weight lists
 - double weights: double data [...]: Weight lists as doubles
 - float weights: float data [...]: Weight lists as floats
```
Weight transpose
--------------------
```
uint64_t offset[chunks-1]: Offset of the file where the transposed weights are
Depending on the type of weights:
 - unsigned/signed weights: varint data [...]: Varint encoded weight lists
 - double weights: double data [...]: Weight lists as doubles
 - float weights: float data [...]: Weight lists as floats
```
