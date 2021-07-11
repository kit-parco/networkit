// no-networkit-format
/*
 * Assortativity.cpp
 *
 *  Created on: Jun 13, 2015
 *      Author: Christian Staudt
 */

#include <cmath>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/correlation/Assortativity.hpp>

namespace NetworKit {

Assortativity::Assortativity(const Graph& G, const std::vector<double>& attribute) : Algorithm(), G(&G), emptyVector(), emptyPartition(), attribute(&attribute), partition(&emptyPartition), nominal(false) {
    if (attribute.size() < G.upperNodeIdBound()) {
        throw std::runtime_error("attribute list has incorrect length: there must be an entry for each node");
    }
}

Assortativity::Assortativity(const Graph& G, const Partition& partition) : Algorithm(), G(&G), emptyVector(), emptyPartition(), attribute(&emptyVector), partition(&partition), nominal(true) {
    if (partition.numberOfElements() < G.upperNodeIdBound()) {
        throw std::runtime_error("partition has incorrect length: there must be an entry for each node");
    }
}


void Assortativity::run() {
    if (nominal) {
        // compact partition so matrix doesn't get unnecessarily large
        Partition P = *partition;
        P.compact();
        // create kxk matrix with entries $e_{ij}$, the fraction of edges connecting nodes of type i to nodes of type j
        count k = P.upperBound();
        std::vector<std::vector<double>> E(k, std::vector<double>(k, 0.0));
        G->forEdges([&](node u, node v) {
            E[P[u]][P[v]] += 1;
        });
        // row and column sums $a_i$ and $b_i$
        std::vector<double> a(k, 0.0);
        std::vector<double> b(k, 0.0);
        // normalize and calculate sums
        count m = G->numberOfEdges();
        for (index i = 0; i < k; ++i) {
            for (index j = 0; j < k; ++j) {
                E[i][j] = E[i][j] / m;
                a[i] += E[i][j];
                b[j] += E[i][j];
            }
        }
        // calculate coefficient $r = \frac{\sum_i e_{ii} - \sum_i a_i b_i}{1 - \sum_i a_i b_i}$
        double diagSum = 0.0;
        double abSum = 0.0;
        for (index i = 0; i < k; ++i) {
            diagSum += E[i][i];
            abSum += a[i] * b[i];
        }
        double r = (diagSum - abSum) / (double) (1 - abSum);
        INFO("diagSum: ", diagSum);
        INFO("abSum: ", abSum);
        INFO("r: ", r);
        coefficient = r;
    } else {
        // assortativity with respect to a continuous, ordinal-scaled node attribute is simply the Pearson correlation coefficient of the lists x and y
        // where (x_u, y_v) are the attributes of connected pairs of nodes
        // r_{xy} := \frac{\sum_{i=1}^n(x_i-\bar x)(y_i-\bar y)}{\sqrt{\sum_{i=1}^n(x_i-\bar x)^2\cdot \sum_{i=1}^n(y_i-\bar y)^2}}

        count m = G->numberOfEdges();
        double xSum = 0.0;
        double ySum = 0.0;
        G->forEdges([&](node u, node v) {
            xSum += (*attribute)[u];
            ySum += (*attribute)[v];
        });
        double xMean = xSum / m;
        double yMean = ySum / m;

        double A = 0.0;
        double B = 0.0;
        double C = 0.0;
        G->forEdges([&](node u, node v) {
            double x = ((*attribute)[u] - xMean);
            double y = ((*attribute)[v] - yMean);
            A +=  x * y;
            B += x * x;
            C += y * y;
        });

        double r = A / sqrt(B * C);
        coefficient = r;
    }
    hasRun = true;
}


double Assortativity::getCoefficient() const {
    return coefficient;
}

std::string Assortativity::toString() const {
    return (std::string)"Assortativity("+((nominal)?"nominal":"ordinal")+")";
}

bool Assortativity::isParallel() const {
    return false;
}


}
