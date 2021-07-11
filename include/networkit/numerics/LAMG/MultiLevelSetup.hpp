// no-networkit-format
/*
 * MultiLevelSetup.hpp
 *
 *  Created on: 10.01.2015
 *      Author: Michael Wegner
 */

#ifndef NETWORKIT_NUMERICS_LAMG_MULTI_LEVEL_SETUP_HPP_
#define NETWORKIT_NUMERICS_LAMG_MULTI_LEVEL_SETUP_HPP_

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/numerics/LAMG/LevelHierarchy.hpp>
#include <networkit/numerics/Smoother.hpp>

#include <limits>
#include <cmath>
#include <vector>

namespace NetworKit {

/**
 * @ingroup
 * Implements the setup phase of LAMG (Lean Algebraic Multigrid by Livne et al.).
 */
template<class Matrix>
class MultiLevelSetup {

#define UNDECIDED std::numeric_limits<index>::max()

private:
    const Smoother<Matrix> &smoother;

    /**
     * Elimination phase of LAMG for the specified Laplacian matrix @a matrix. If an elemination stage is created,
     * it is stored in the LevelHierarchy @a hierarchy and the method returns @code{true}. Otherwise, @code{false}
     * is returned.
     * @param matrix Laplacian matrix for which an elimination stage should be created.
     * @param hierarchy LevelHierarchy in which the created elimination stage is stored (if created)
     * @return @code{True} if an elimination stage has been created, otherwise@code{false}.
     */
    bool coarseningElimination(Matrix& matrix, LevelHierarchy<Matrix>& hierarchy) const;

    /**
     * Scans the Laplacian matrix for nodes with a low degree (i.e. nodes with less than 5 neighbors). For each node
     * having low degree and independent from already found low degree nodes, @code{true} is stored in @a fNode. The
     * @a stage parameter specifies if we are in the first or subsequent stages during elimination.
     * @param matrix Laplacian matrix.
     * @param fNode[out] For each node, @code{true} if the node is of low degree and @code{false} otherwise.
     * @param stage The stage of the elimination phase.
     * @return Number of low degree nodes found.
     */
    count lowDegreeSweep(const Matrix& matrix, std::vector<bool>& fNode, index stage) const;

    /**
     * Computes the projection matrix @a P and the @a q vector used to restrict and interpolate the matrix for an
     * elimination stage.
     * @param matrix Laplacian matrix.
     * @param fSet Set of nodes having low degree.
     * @param coarseIndex Set of nodes equal to V \setminus fSet
     * @param P[out] The projection matrix.
     * @param q[out] The q vector.
     */
    void eliminationOperators(const Matrix& matrix, const std::vector<index>& fSet, const std::vector<index>& coarseIndex, Matrix& P, Vector& q) const;

    /**
     * Aggregation phase of LAMG for the specified Laplacian matrix @a matrix. The coarsened matrix is stored in the
     * LevelHierarchy @a hierarchy. The test vector @a tv is used for determining the affinity of nodes. The parameter
     * @a numTVVectors specifies how many test vectors will be created to determine the affinities between nodes.
     * @param matrix Laplacian matrix.
     * @param hierarchy LevelHierarchy to store the stage.
     * @param tv Test vector.
     * @param numTVVectors Number of test vectors to use.
     */
    void coarseningAggregation(Matrix& matrix, LevelHierarchy<Matrix>& hierarchy, Vector& tv, count numTVVectors) const;

    /**
     * Create @a numVectors test vectors for the given Laplacian matrix @a matrix. The test vector @a tv will be
     * reused.
     * @param matrix Laplacian matrix.
     * @param tv Test vector.
     * @param numVectors Number of test vectors to create.
     * @return The created test vectors.
     */
    std::vector<Vector> generateTVs(const Matrix &matrix, Vector &tv, count numVectors) const;

    /**
     * Adds high degree nodes as seeds to @a status.
     * @param matrix Laplacian matrix.
     * @param status[out] High degree nodes are flagged as seed.
     */
    void addHighDegreeSeedNodes(const Matrix& matrix, std::vector<index>& status) const;

    /**
     * Aggregates loose nodes (nodes with low adjacency) together.
     * @param strongAdjMatrix Strong adjacency matrix.
     * @param status[out] Aggregates loose nodes together and labels them in @a status accordingly.
     * @param nc[out] The altered number of coarse nodes.
     */
    void aggregateLooseNodes(const Matrix& strongAdjMatrix, std::vector<index>& status, count& nc) const;

    /**
     * Computes the strong adjacency matrix for the given Laplacian matrix @a matrix.
     * @param matrix Laplacian matrix.
     * @param strongAdjMatrix[out] The resulting strong adjacency matrix.
     */
    void computeStrongAdjacencyMatrix(const Matrix& matrix, Matrix& strongAdjMatrix) const;

    /**
     * Computes the affinity matrix for the given Laplacian matrix @a matrix and the test vectors @a tVs.
     * @param matrix Laplacian matrix.
     * @param tVs Test vectors.
     * @param affinityMatrix[out] The resulting affinity matrix.
     */
    void computeAffinityMatrix(const Matrix& matrix, const std::vector<Vector>& tVs, Matrix& affinityMatrix) const;

    /**
     * Models one stage in the aggregation phase. New aggregates are labeled accordingly in @a status.
     * @param matrix Laplacian matrix.
     * @param nc Number of coarse nodes.
     * @param strongAdjMatrix Strong adjacency matrix.
     * @param affinityMatrix Affinity matrix.
     * @param tVs[out] Test vectors.
     * @param status[out] Aggregation labels.
     */
    void aggregationStage(const Matrix& matrix, count& nc, const Matrix& strongAdjMatrix, const Matrix& affinityMatrix, std::vector<Vector>& tVs, std::vector<index>& status) const;

    /**
     * Computes strong (cf. LAMG paper) neighbors and stores them sorted into @a bins.
     * @param affinityMatrix Affinity matrix.
     * @param status Aggregation labels.
     * @param bins[out] Strong neighbors sorted into bins.
     */
    void computeStrongNeighbors(const Matrix& affinityMatrix, const std::vector<index>& status, std::vector<std::vector<index>>& bins) const;

    /**
     * Finds the best seed for node @a u and stores it in @a s.
     * @param strongAdjMatrix
     * @param affinityMatrix Affinity matrix.
     * @param diag Vector of diagonal entries of the Laplacian matrix.
     * @param tVs Test vectors.
     * @param status Aggregation labels.
     * @param u The node to find the best seed for.
     * @param s[out] The best seed for node @a u.
     * @return @code{True} if a seed has been found for @a u, @code{false} otherwise.
     */
    bool findBestSeedEnergyCorrected(const Matrix &strongAdjMatrix, const Matrix &affinityMatrix,
                                     const std::vector<double> &diag,
                                     const std::vector<Vector> &tVs,
                                     const std::vector<index> &status, index u, index &s) const;

    /**
     * Determines if the Laplacian matrix @a A can be coarsened further.
     * @param A Laplacian matrix.
     * @return @code{True} if @a A can be coarsened further, @code{false} otherwise.
     */
    inline bool canCoarsen(const Matrix& A) const {
        return A.numberOfRows() > MAX_DIRECT_SOLVE_SIZE;
    }

    /**
     * Determines if the relaxation is fast enough to stop coarsening.
     * @param A Laplacian matrix.
     * @param lvlIndex The number of levels already created in the hierarchy.
     * @param tv Test vector.
     * @return @code{True} if convergence of relaxation is fast, @code{false} otherwise.
     */
    bool isRelaxationFast(const Matrix& A, index lvlIndex, Vector& tv) const;

    /**
     * Computes the coarsened matrix of @a matrix by means of the projection matrix @a P and stores the result in @a B.
     * @param P Projection matrix.
     * @param A Laplacian matrix.
     * @param PColIndex Stores the column index of the 1 entry at each row.
     * @param PRowIndex Stores the row index of the 1 entry at each column.
     * @param B[out] Resulting coarsened Laplacian matrix.
     */
    void galerkinOperator(const Matrix& P, const Matrix& A, const std::vector<index>& PColIndex, const std::vector<std::vector<index>>& PRowIndex, Matrix& B) const;

    /**
     * Creates a @a hierarchy for the given Laplacian matrix @a matrix.
     * @param matrix Laplcian matrix.
     * @param hierarchy[out] The constructed hierarchy.
     */
    void setupForMatrix(Matrix& matrix, LevelHierarchy<Matrix>& hierarchy) const;

public:
    /**
     * Creates an instance of MultiLevelSetup with the specified @a smoother used for relaxing during the setup phase.
     * @param smoother Reference to smoother.
     */
    MultiLevelSetup(const Smoother<Matrix>& smoother) : smoother(smoother) {}

    /**
     * Creates a @ hierarchy for the given Laplacian matrix of the graph @a G.
     * @param G The graph.
     * @param hierarchy[out] The constructed hierarchy.
     */
    void setup(const Graph& G, LevelHierarchy<Matrix>& hierarchy) const {
        setup(Matrix::laplacianMatrix(G), hierarchy);
    }

    /**
     * Creates a @a hierarchy for the given Laplacian matrix @a matrix.
     * @param matrix Laplcian matrix.
     * @param hierarchy[out] The constructed hierarchy.
     */
    void setup(const Matrix& matrix, LevelHierarchy<Matrix>& hierarchy) const;
};

template<class Matrix>
void MultiLevelSetup<Matrix>::setup(const Matrix& matrix, LevelHierarchy<Matrix>& hierarchy) const {
    CSRMatrix A = matrix;
    setupForMatrix(A, hierarchy);
}

template<class Matrix>
void MultiLevelSetup<Matrix>::setupForMatrix(Matrix& A, LevelHierarchy<Matrix>& hierarchy) const {
    hierarchy.addFinestLevel(A);

#ifndef NDEBUG
    DEBUG("FINEST\t", A.numberOfRows(), "\t", A.nnz());
#endif // NDEBUG

    bool doneCoarsening = false;
    count numTVs = TV_NUM;
    index level = 0;

    while (!doneCoarsening) {
        // ELIMINATION
        if (coarseningElimination(A, hierarchy)) {
            if (!canCoarsen(A)) doneCoarsening = true;
            level++;
#ifndef NDEBUG
            DEBUG(level, " ELIM\t\t", A.numberOfRows(), "\t", A.nnz() / 2);
#endif // NDEBUG
        }

        // AGGREGATION
        Vector tv;
        if (doneCoarsening || isRelaxationFast(A, level, tv)) {
            doneCoarsening = true;
        } else {
            coarseningAggregation(A, hierarchy, tv, numTVs);
            level++;
#ifndef NDEBUG
            DEBUG(level, " AGG\t\t", A.numberOfRows(), "\t", A.nnz() / 2);
#endif // NDEBUG
            if (numTVs < TV_MAX) {
                numTVs += TV_INC;
            }
        }

        if (!canCoarsen(A)) doneCoarsening = true;
    }

    hierarchy.setLastAsCoarsest();
}


template<class Matrix>
bool MultiLevelSetup<Matrix>::coarseningElimination(Matrix& matrix, LevelHierarchy<Matrix>& hierarchy) const {
    std::vector<EliminationStage<Matrix>> coarseningStages;
    count stageNum = 0;
    while (stageNum < SETUP_ELIMINATION_MAX_STAGES) {
        if (matrix.numberOfRows() <= MAX_DIRECT_SOLVE_SIZE) break; // we do not need to coarsen the matrix any further
        std::vector<bool> fNode;
        count nf = lowDegreeSweep(matrix, fNode, stageNum);
        count nc = matrix.numberOfRows() - nf;

        if (nc == 0) { // do not eliminate all nodes -> leave one entry in c
            nc = 1;
            nf--;
        }

        // add f nodes to fSet and c nodes to cSet
        std::vector<index> fSet(nf);
        std::vector<index> cSet(nc);

        std::vector<index> coarseIndex(matrix.numberOfRows());
        count numFNodes = 0;
        for (index i = 0, fIndex = 0, cIndex = 0; i < matrix.numberOfRows(); ++i) {
            if (fNode[i] && fIndex < nf) {
                coarseIndex[i] = fIndex;
                fSet[fIndex++] = i;
                numFNodes++;
            } else {
                coarseIndex[i] = cIndex;
                cSet[cIndex++] = i;
            }
        }

        if (nf <= SETUP_ELIMINATION_MIN_ELIM_FRACTION * matrix.numberOfRows()) {
            break;
        }

        Matrix P;
        Vector q;
        eliminationOperators(matrix, fSet, coarseIndex, P, q);
        coarseningStages.push_back(EliminationStage<Matrix>(P, q, fSet, cSet));

        Matrix Acc = matrix.extract(cSet, cSet); // Schur complement
        Matrix Acf = matrix.extract(cSet, fSet); // Schur complement

        matrix = Acc + Acf * P;
        stageNum++;
    }

    if (stageNum != 0) { // we have coarsened the matrix
        hierarchy.addEliminationLevel(matrix, coarseningStages);
        return true;
    }

    return false;
}

template<class Matrix>
count MultiLevelSetup<Matrix>::lowDegreeSweep(const Matrix& matrix, std::vector<bool>& fNode, index stage) const {
    fNode.resize(matrix.numberOfRows(), true); // first mark all nodes as f nodes
    count numFNodes = 0;
    int degreeOffset = stage != 0;

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        if ((int) matrix.nnzInRow(i) - degreeOffset <= (int)SETUP_ELIMINATION_MAX_DEGREE && fNode[i]) { // node i has degree <= 4 and can be eliminated
            numFNodes++;
            matrix.forNonZeroElementsInRow(i, [&](index j, edgeweight /*w*/){ // to maintain independence, mark all neighbors as not eliminated
                if (j != i)	{ // all neighbors of this f node are c nodes
                    fNode[j] = false;
                }
            });
        } else { // node has high degree, thus it is a c node
            fNode[i] = false;
        }
    }

    return numFNodes;
}

template<class Matrix>
void MultiLevelSetup<Matrix>::eliminationOperators(const Matrix& matrix, const std::vector<index>& fSet, const std::vector<index>& coarseIndex, Matrix& P, Vector& q) const {
    std::vector<Triplet> triples;
    q = Vector(fSet.size());
    for (index k = 0; k < fSet.size(); ++k) { // Afc
        matrix.forNonZeroElementsInRow(fSet[k], [&](index j, edgeweight w){
            if (fSet[k] == j) {
                q[k] = 1.0 / w;
            } else {
                triples.push_back({k, coarseIndex[j], w});
            }
        });
    }

    for (index i = 0; i < triples.size(); ++i) { // * -Aff^-1
        triples[i].value *= -q[triples[i].row];
    }

    P = Matrix(fSet.size(), coarseIndex.size() - fSet.size(), triples);
}



template<class Matrix>
void MultiLevelSetup<Matrix>::coarseningAggregation(Matrix& matrix, LevelHierarchy<Matrix>& hierarchy, Vector& tv, count numTVVectors) const {
    Vector B(SETUP_MAX_AGGREGATION_STAGES, std::numeric_limits<double>::max());
    std::vector<std::vector<index>> S(SETUP_MAX_AGGREGATION_STAGES, std::vector<index>(matrix.numberOfRows(), std::numeric_limits<index>::max()));
    std::vector<index> status(matrix.numberOfRows(), UNDECIDED);
    std::vector<count> nc(SETUP_MAX_AGGREGATION_STAGES, matrix.numberOfRows());

    double alpha = 1.0;
    double maxCoarseningRatio = SETUP_COARSENING_WORK_GUARD / SETUP_CYCLE_INDEX;
    count stage = 0;
    count nC = matrix.numberOfRows();


    // generate TVs
    std::vector<Vector> tVs = generateTVs(matrix, tv, numTVVectors);

    // compute strong adjacency matrix
    Matrix Wstrong;
    computeStrongAdjacencyMatrix(matrix, Wstrong);

    // compute affinityMatrix
    Matrix affinityMatrix;
    computeAffinityMatrix(Wstrong, tVs, affinityMatrix);

    // mark all locally high-degree nodes as seeds
    addHighDegreeSeedNodes(matrix, status);

    // aggregate all loose nodes
    aggregateLooseNodes(Wstrong, status, nC);

    nc[0] = nC;
    while (stage < SETUP_MIN_AGGREGATION_STAGES || (alpha >= maxCoarseningRatio && stage < SETUP_MAX_AGGREGATION_STAGES)) {
        nC = stage > 0? nc[stage - 1] : nc[0];

        // aggregation stage
        aggregationStage(matrix, nC, Wstrong, affinityMatrix, tVs, status);

        alpha = (double) nC / (double) matrix.numberOfRows();
        alpha <= maxCoarseningRatio? B[stage] = 1.0-alpha : B[stage] = 1.0+alpha;

        S[stage] = status;
        nc[stage] = nC;
        stage++;
    }

    double min = B[0];
    index bestAggregate = 0;
    for (index i = 1; i < stage; ++i) {
        if (B[i] < min) {
            bestAggregate = i;
            min = B[i];
        }
    }

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        if (S[bestAggregate][i] == UNDECIDED) { // undediced nodes become their own seeds
            S[bestAggregate][i] = i;
        }
    }

    std::vector<index> indexFine(matrix.numberOfRows(), 0);
    index newIndex = 0;
    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        if (S[bestAggregate][i] == i) {
            indexFine[i] = newIndex++;
        }
    }

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        status[i] = indexFine[S[bestAggregate][i]];
    }

    assert(newIndex == nc[bestAggregate]);

    // create interpolation matrix
    std::vector<Triplet> pTriples(matrix.numberOfRows());
    std::vector<Triplet> rTriples(matrix.numberOfRows());
    std::vector<index> PColIndex(matrix.numberOfRows());
    std::vector<std::vector<index>> PRowIndex(nc[bestAggregate]);

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        pTriples[i] = {i, status[i], 1};
        rTriples[i] = {status[i], i, 1};
        PColIndex[i] = status[i];
        PRowIndex[status[i]].push_back(i);
    }

    Matrix P(matrix.numberOfRows(), nc[bestAggregate], pTriples);
    Matrix R(nc[bestAggregate], matrix.numberOfRows(), rTriples);

    // create coarsened laplacian
    galerkinOperator(P, matrix, PColIndex, PRowIndex, matrix);

    hierarchy.addAggregationLevel(matrix, P, R);
}

template<class Matrix>
std::vector<Vector> MultiLevelSetup<Matrix>::generateTVs(const Matrix& matrix, Vector& tv, count numVectors) const {
    std::vector<Vector> testVectors(numVectors, Vector(matrix.numberOfColumns()));

    testVectors[0] = tv;

    if (numVectors > 1) {
        Vector b(matrix.numberOfColumns(), 0.0);
#pragma omp parallel for
        for (omp_index i = 1; i < static_cast<omp_index>(numVectors); ++i) {
            for (count j = 0; j < matrix.numberOfColumns(); ++j) {
                testVectors[i][j] = 2 * Aux::Random::probability() - 1;
            }

            testVectors[i] = smoother.relax(matrix, b, testVectors[i], SETUP_TV_SWEEPS);
        }
    }

    return testVectors;
}

template<class Matrix>
void MultiLevelSetup<Matrix>::addHighDegreeSeedNodes(const Matrix& matrix, std::vector<index>& status) const {
    std::vector<count> deg(matrix.numberOfRows());
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        deg[i] = matrix.nnzInRow(i) - 1;
    }

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        double num = 0.0;
        double denom = 0.0;
        matrix.forNonZeroElementsInRow(i, [&](index j, double value){
            if (i != j) {
                num += std::abs(value) * (double) deg[j];
            } else {
                denom = std::abs(value);
            }
        });


        if ((double) deg[i] >= SETUP_AGGREGATION_DEGREE_THRESHOLD * (num / denom)) { // high degree node becomes seed
            status[i] = i;
        }
    }
}

template<class Matrix>
void MultiLevelSetup<Matrix>::aggregateLooseNodes(const Matrix& strongAdjMatrix, std::vector<index>& status, count& nc) const {
    std::vector<index> looseNodes;
    for (index i = 0; i < strongAdjMatrix.numberOfRows(); ++i) {
        double max = std::numeric_limits<double>::min();
        strongAdjMatrix.forNonZeroElementsInRow(i, [&](index /*j*/, double value) {
            if (value > max) max = value;
        });

        if (std::abs(max) < 1e-9 || max == std::numeric_limits<double>::min()) {
            looseNodes.push_back(i);
        }
    }

    if (!looseNodes.empty()) {
        status[looseNodes[0]] = looseNodes[0]; // mark first as seed
        for (index k = 1; k < looseNodes.size(); ++k) {
            status[looseNodes[k]] = looseNodes[0]; // first loose nodes becomes seed
        }

        nc -= looseNodes.size() - 1;
    }
}



template<class Matrix>
void MultiLevelSetup<Matrix>::computeStrongAdjacencyMatrix(const Matrix& matrix, Matrix& strongAdjMatrix) const {
    std::vector<double> maxNeighbor(matrix.numberOfRows(), std::numeric_limits<double>::min());
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        matrix.forNonZeroElementsInRow(i, [&](index j, double value) {
            if (i != j && -value > maxNeighbor[i]) {
                maxNeighbor[i] = -value;
            }
        });
    }

    std::vector<index> rowIdx(matrix.numberOfRows()+1, 0);

    matrix.parallelForNonZeroElementsInRowOrder([&](index i, index j, double value) {
        if (i != j && std::abs(value) >= 0.1 * std::min(maxNeighbor[i], maxNeighbor[j])) {
            ++rowIdx[i+1];
        }
    });

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        rowIdx[i+1] += rowIdx[i];
    }

    count nnz = rowIdx[matrix.numberOfRows()];
    std::vector<Triplet> triplets(nnz);

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        index cIdx = rowIdx[i];
        matrix.forNonZeroElementsInRow(i, [&](index j, double value) {
            if (i != j && std::abs(value) >= 0.1 * std::min(maxNeighbor[i], maxNeighbor[j])) {
                triplets[cIdx] = {static_cast<index>(i),j,-value};
                ++cIdx;
            }
        });
    }

    strongAdjMatrix = Matrix(matrix.numberOfRows(), matrix.numberOfColumns(), triplets);
}


template<class Matrix>
void MultiLevelSetup<Matrix>::computeAffinityMatrix(const Matrix& matrix, const std::vector<Vector>& tVs, Matrix& affinityMatrix) const {
    assert(!tVs.empty());

    std::vector<index> rowIdx(matrix.numberOfRows()+1);
    std::vector<Triplet> triplets(matrix.nnz());

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        rowIdx[i+1] = matrix.nnzInRow(i);
    }

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        rowIdx[i+1] += rowIdx[i];
    }

    std::vector<double> normSquared(matrix.numberOfRows(), 0.0);
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        for (index k = 0; k < tVs.size(); ++k) {
            normSquared[i] += tVs[k][i] * tVs[k][i];
        }
    }

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        double nir = 1.0 / normSquared[i];
        index cIdx = rowIdx[i];
        matrix.forNonZeroElementsInRow(i, [&](index j, double /*val*/) {
            double ij = 0.0;
            for (index k = 0; k < tVs.size(); ++k) {
                ij += tVs[k][i] * tVs[k][j];
            }

            double value = (ij * ij) * nir / normSquared[j];
            triplets[cIdx] = {static_cast<index>(i),j,value};
            ++cIdx;
        });
    }

    affinityMatrix = Matrix(matrix.numberOfRows(), matrix.numberOfColumns(), triplets);
}

template<class Matrix>
void MultiLevelSetup<Matrix>::aggregationStage(const Matrix& matrix, count& nc, const Matrix& strongAdjMatrix, const Matrix& affinityMatrix, std::vector<Vector>& tVs, std::vector<index>& status) const {
    std::vector<std::vector<index>> bins(10);
    computeStrongNeighbors(affinityMatrix, status, bins);

    std::vector<double> diag(matrix.numberOfRows(), 0.0);
#pragma omp parallel for
    for (omp_index i = 0 ; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        diag[i] = matrix(i,i);
    }

    for (index k = bins.size(); k-- > 0;) { // iterate over undecided nodes with strong neighbors in decreasing order of strongest neighbor
        for (index i : bins[k]) {
            if (status[i] == UNDECIDED) { // node is still undecided
                index s = 0;
                if (findBestSeedEnergyCorrected(strongAdjMatrix, affinityMatrix, diag, tVs, status, i, s)) {
                    status[s] = s; // s becomes seed
                    status[i] = s; // i's seed is s
                    nc--;

                    for (index j = 0; j < tVs.size(); ++j) { // update test vectors
                        tVs[j][i] = tVs[j][s];
                    }
                }
            }
        }

        if (nc <= matrix.numberOfRows() * SETUP_COARSENING_WORK_GUARD / SETUP_CYCLE_INDEX) {
            break;
        }
    } // iterate over bins
}

template<class Matrix>
void MultiLevelSetup<Matrix>::computeStrongNeighbors(const Matrix& affinityMatrix, const std::vector<index>& status, std::vector<std::vector<index>>& bins) const {
    std::vector<bool> undecided(affinityMatrix.numberOfRows(), false);
    std::vector<double> maxNeighbor(affinityMatrix.numberOfRows(), std::numeric_limits<double>::min());
    double overallMax = 0.0;
    double overallMin = std::numeric_limits<double>::max();

    affinityMatrix.parallelForNonZeroElementsInRowOrder([&](index i, index j, double value) { // determine the highest affinity neighbor of each node
        if (status[i] == UNDECIDED && (status[j] == UNDECIDED || status[j] == j)) { // i is UNDECIDED and its neighbor j is also UNDECIDED or SEED
            if (value > maxNeighbor[i]) {
                maxNeighbor[i] = value;
            }
            undecided[i] = true;
        }
    });

    for (index i = 0; i < affinityMatrix.numberOfRows(); ++i) {
        if (maxNeighbor[i] > overallMax) {
            overallMax = maxNeighbor[i];
        }
        if (maxNeighbor[i] < overallMin) {
            overallMin = maxNeighbor[i];
        }
    }

    double h = fabs(overallMax - overallMin) < 1e-15? 1.0 : (double) bins.size() / (overallMax - overallMin);
    for (index i = 0; i < affinityMatrix.numberOfRows(); ++i) {
        if (undecided[i]) { // undecided nodes with strong neighbors
            index binIndex = (index) std::floor(h * (maxNeighbor[i] - overallMin));
            if (binIndex == bins.size()) { // last interval is closed on the right
                binIndex--;
            }

            assert(binIndex < bins.size());
            bins[binIndex].push_back(i);
        }
    }
}

template<class Matrix>
bool MultiLevelSetup<Matrix>::findBestSeedEnergyCorrected(const Matrix& strongAdjMatrix, const Matrix& affinityMatrix, const std::vector<double>& diag, const std::vector<Vector>& tVs, const std::vector<index>& status, const index u, index& s) const {
    bool foundSeed = false;
    std::vector<double> r(tVs.size(), 0.0);
    std::vector<double> q(tVs.size(), 0.0);
    std::vector<double> E(tVs.size(), 0.0);

    double d = diag[u];
    double d2 = 0.5 * diag[u];
    for (index k = 0; k < tVs.size(); ++k) {
        double rr = 0.0;
        double qq = 0.0;
        strongAdjMatrix.forNonZeroElementsInRow(u, [&](index v, double value) {
            rr += value * tVs[k][v];
            qq += value * 0.5 * tVs[k][v] * tVs[k][v];
        });

        r[k] = rr;
        q[k] = qq;
        double y = rr/d;
        E[k] = (d2*y - rr)*y + qq;
    }

    double maxNeighbor = -1.0;
    affinityMatrix.forNonZeroElementsInRow(u, [&](index v, double value) {
        if (status[v] == UNDECIDED || status[v] == v) {
            double maxMu = -1.0;
            bool smallRatio = true;
            for (index k = 0; k < tVs.size(); ++k) {
                double xv = tVs[k][v];
                double Ec = (d2*xv - r[k])*xv + q[k];
                double mu = Ec / (E[k] + 1e-15);

                if (mu > maxMu) {
                    maxMu = mu;
                }
                if (maxMu > 2.5) {
                    smallRatio = false;
                    break;
                }
            }

            if (smallRatio && value > maxNeighbor) {
                maxNeighbor = value;
                s = v;
                foundSeed = true;
            }
        }
    });

    return foundSeed;
}

template<class Matrix>
bool MultiLevelSetup<Matrix>::isRelaxationFast(const Matrix& A, index lvlIndex, Vector& tv) const {
    count nu = SETUP_RELAX_ACF_MIN_SWEEPS + 2 * (lvlIndex - 1);
    count tvNu = SETUP_TV_SWEEPS;
    count initial = 3;

    // create testVector in [-1,1]
    tv = Vector(A.numberOfRows());
    for (index i = 0; i < tv.getDimension(); ++i) {
        tv[i] = 2.0 * Aux::Random::probability() - 1.0;
    }

    Vector b(A.numberOfRows(), 0.0);
    Vector x = tv;
    x = smoother.relax(A, b, x, initial);
    tv = smoother.relax(A, b, x, tvNu - initial);
    Vector y = smoother.relax(A, b, tv, nu - tvNu);
    double relaxAcf = std::pow((y - y.mean()).length() / (x - x.mean()).length(), (double) 1.0 / (double) (nu - initial));
    return relaxAcf <= SETUP_MAX_COARSE_RELAX_ACF || !canCoarsen(A);
}



template<class Matrix>
void MultiLevelSetup<Matrix>::galerkinOperator(const Matrix& P, const Matrix& A, const std::vector<index>& PColIndex, const std::vector<std::vector<index>>& PRowIndex, Matrix& B) const {
    std::vector<Triplet> triplets;
    SparseAccumulator spa(P.numberOfColumns());
    for (index i = 0; i < P.numberOfColumns(); ++i) {
        for (index k : PRowIndex[i]) {
            double Pki = P(k,i);
            A.forNonZeroElementsInRow(k, [&](index l, double value) {
                index j = PColIndex[l];
                spa.scatter(Pki * value * P(l, j), j);
            });
        }

        spa.gather([&](index i, index j, double value) {
            triplets.push_back({i,j,value});
        });

        spa.increaseRow();
    }

    B = CSRMatrix(P.numberOfColumns(), P.numberOfColumns(), triplets);
}

template<>
inline void MultiLevelSetup<CSRMatrix>::setup(const CSRMatrix& matrix, LevelHierarchy<CSRMatrix>& hierarchy) const {
    CSRMatrix A = matrix;
    A.sort();
    setupForMatrix(A, hierarchy);
}

template<>
inline void MultiLevelSetup<CSRMatrix>::computeStrongAdjacencyMatrix(const CSRMatrix& matrix, CSRMatrix& strongAdjMatrix) const {
    std::vector<double> maxNeighbor(matrix.numberOfRows(), std::numeric_limits<double>::min());
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        matrix.forNonZeroElementsInRow(i, [&](index j, double value) {
            if (i != j && -value > maxNeighbor[i]) {
                maxNeighbor[i] = -value;
            }
        });
    }

    std::vector<index> rowIdx(matrix.numberOfRows()+1, 0);
    matrix.parallelForNonZeroElementsInRowOrder([&](index i, index j, double value) {
        if (i != j && std::abs(value) >= 0.1 * std::min(maxNeighbor[i], maxNeighbor[j])) {
            ++rowIdx[i+1];
        }
    });

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        rowIdx[i+1] += rowIdx[i];
    }

    count nnz = rowIdx[matrix.numberOfRows()];
    std::vector<index> columnIdx(nnz);
    std::vector<double> nonZeros(nnz);

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        index cIdx = rowIdx[i];
        matrix.forNonZeroElementsInRow(i, [&](index j, double value) {
            if (i != j && std::abs(value) >= 0.1 * std::min(maxNeighbor[i], maxNeighbor[j])) {
                columnIdx[cIdx] = j;
                nonZeros[cIdx] = -value;
                ++cIdx;
            }
        });
    }

    strongAdjMatrix = CSRMatrix(matrix.numberOfRows(), matrix.numberOfColumns(), rowIdx, columnIdx, nonZeros, 0.0, matrix.sorted());
}

template<>
inline void MultiLevelSetup<CSRMatrix>::computeAffinityMatrix(const CSRMatrix& matrix, const std::vector<Vector> &tVs, CSRMatrix& affinityMatrix) const {
    assert(!tVs.empty());

    std::vector<index> rowIdx(matrix.numberOfRows()+1);
    std::vector<index> columnIdx(matrix.nnz());
    std::vector<double> nonZeros(matrix.nnz());

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        rowIdx[i+1] = matrix.nnzInRow(i);
    }

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        rowIdx[i+1] += rowIdx[i];
    }

    std::vector<double> normSquared(matrix.numberOfRows(), 0.0);
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        for (index k = 0; k < tVs.size(); ++k) {
            normSquared[i] += tVs[k][i] * tVs[k][i];
        }
    }

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
        double nir = 1.0 / normSquared[i];
        index cIdx = rowIdx[i];
        matrix.forNonZeroElementsInRow(i, [&](index j, double /*val*/) {
            double ij = 0.0;
            for (index k = 0; k < tVs.size(); ++k) {
                ij += tVs[k][i] * tVs[k][j];
            }

            double value = (ij * ij) * nir / normSquared[j];
            columnIdx[cIdx] = j;
            nonZeros[cIdx] = value;
            ++cIdx;
        });
    }

    affinityMatrix = CSRMatrix(matrix.numberOfRows(), matrix.numberOfColumns(), rowIdx, columnIdx, nonZeros, 0.0, matrix.sorted());
}

template<>
inline void MultiLevelSetup<CSRMatrix>::galerkinOperator(const CSRMatrix& P, const CSRMatrix& A, const std::vector<index>& PColIndex, const std::vector<std::vector<index>>& PRowIndex, CSRMatrix& B) const {
    std::vector<Triplet> triplets;
    SparseAccumulator spa(P.numberOfColumns());
    for (index i = 0; i < P.numberOfColumns(); ++i) {
        for (index k : PRowIndex[i]) {
            double Pki = P(k,i);
            A.forNonZeroElementsInRow(k, [&](index l, double value) {
                index j = PColIndex[l];
                spa.scatter(Pki * value * P(l, j), j);
            });
        }

        spa.gather([&](index i, index j, double value) {
            triplets.push_back({i,j,value});
        });

        spa.increaseRow();
    }

    B = CSRMatrix(P.numberOfColumns(), P.numberOfColumns(), triplets, 0.0, true);
}


template<>
inline void MultiLevelSetup<CSRMatrix>::eliminationOperators(const CSRMatrix& matrix, const std::vector<index>& fSet, const std::vector<index>& coarseIndex, CSRMatrix& P, Vector& q) const {
    std::vector<Triplet> triples;
    q = Vector(fSet.size());
    for (index k = 0; k < fSet.size(); ++k) { // Afc
        matrix.forNonZeroElementsInRow(fSet[k], [&](index j, edgeweight w){
            if (fSet[k] == j) {
                q[k] = 1.0 / w;
            } else {
                triples.push_back({k, coarseIndex[j], w});
            }
        });
    }

    for (index i = 0; i < triples.size(); ++i) { // * -Aff^-1
        triples[i].value *= -q[triples[i].row];
    }

    P = CSRMatrix(fSet.size(), coarseIndex.size() - fSet.size(), triples, 0.0, matrix.sorted());
}

template<>
inline void MultiLevelSetup<CSRMatrix>::coarseningAggregation(CSRMatrix& matrix, LevelHierarchy<CSRMatrix>& hierarchy, Vector& tv, count numTVVectors) const {
    Vector B(SETUP_MAX_AGGREGATION_STAGES, std::numeric_limits<double>::max());
    std::vector<std::vector<index>> S(SETUP_MAX_AGGREGATION_STAGES, std::vector<index>(matrix.numberOfRows(), std::numeric_limits<index>::max()));
    std::vector<index> status(matrix.numberOfRows(), UNDECIDED);
    std::vector<count> nc(SETUP_MAX_AGGREGATION_STAGES, matrix.numberOfRows());

    double alpha = 1.0;
    double maxCoarseningRatio = SETUP_COARSENING_WORK_GUARD / SETUP_CYCLE_INDEX;
    count stage = 0;
    count nC = matrix.numberOfRows();


    // generate TVs
    std::vector<Vector> tVs = generateTVs(matrix, tv, numTVVectors);

    // compute strong adjacency matrix
    CSRMatrix Wstrong;
    computeStrongAdjacencyMatrix(matrix, Wstrong);

    // compute affinityMatrix
    CSRMatrix affinityMatrix;
    computeAffinityMatrix(Wstrong, tVs, affinityMatrix);

    // mark all locally high-degree nodes as seeds
    addHighDegreeSeedNodes(matrix, status);

    // aggregate all loose nodes
    aggregateLooseNodes(Wstrong, status, nC);

    nc[0] = nC;
    while (stage < SETUP_MIN_AGGREGATION_STAGES || (alpha >= maxCoarseningRatio && stage < SETUP_MAX_AGGREGATION_STAGES)) {
        nC = stage > 0? nc[stage - 1] : nc[0];

        // aggregation stage
        aggregationStage(matrix, nC, Wstrong, affinityMatrix, tVs, status);

        alpha = (double) nC / (double) matrix.numberOfRows();
        alpha <= maxCoarseningRatio? B[stage] = 1.0-alpha : B[stage] = 1.0+alpha;

        S[stage] = status;
        nc[stage] = nC;
        stage++;
    }

    double min = B[0];
    index bestAggregate = 0;
    for (index i = 1; i < stage; ++i) {
        if (B[i] < min) {
            bestAggregate = i;
            min = B[i];
        }
    }

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        if (S[bestAggregate][i] == UNDECIDED) { // undediced nodes become their own seeds
            S[bestAggregate][i] = i;
        }
    }

    std::vector<index> indexFine(matrix.numberOfRows(), 0);
    index newIndex = 0;
    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        if (S[bestAggregate][i] == i) {
            indexFine[i] = newIndex++;
        }
    }

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        status[i] = indexFine[S[bestAggregate][i]];
    }

    assert(newIndex == nc[bestAggregate]);

    // create interpolation matrix
    std::vector<Triplet> pTriples(matrix.numberOfRows());
    std::vector<Triplet> rTriples(matrix.numberOfRows());
    std::vector<index> PColIndex(matrix.numberOfRows());
    std::vector<std::vector<index>> PRowIndex(nc[bestAggregate]);

    for (index i = 0; i < matrix.numberOfRows(); ++i) {
        pTriples[i] = {i, status[i], 1};
        rTriples[i] = {status[i], i, 1};
        PColIndex[i] = status[i];
        PRowIndex[status[i]].push_back(i);
    }

    CSRMatrix P(matrix.numberOfRows(), nc[bestAggregate], pTriples, 0.0, matrix.sorted());
    CSRMatrix R(nc[bestAggregate], matrix.numberOfRows(), rTriples, 0.0, matrix.sorted());

    // create coarsened laplacian
    galerkinOperator(P, matrix, PColIndex, PRowIndex, matrix);

    hierarchy.addAggregationLevel(matrix, P, R);
}

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_LAMG_MULTI_LEVEL_SETUP_HPP_
