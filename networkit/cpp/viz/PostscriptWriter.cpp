/*
 * PostscriptWriter.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/viz/PostscriptWriter.hpp>

namespace NetworKit {

namespace PostscriptWriterColors {
struct RGBColor {
    float r;
    float g;
    float b;
};

static RGBColor fromCyclicRotation(size_t index) {
    constexpr RGBColor colors[] = {
        RGBColor({1.0f, 0.0f, 0.0f}), RGBColor({1.0f, 0.5f, 0.0f}), RGBColor({1.0f, 1.0f, 0.0f}),
        RGBColor({0.5f, 1.0f, 0.0f}), RGBColor({0.0f, 1.0f, 0.0f}), RGBColor({0.0f, 1.0f, 0.5f}),
        RGBColor({0.0f, 1.0f, 1.0f}), RGBColor({0.0f, 0.5f, 1.0f}), RGBColor({0.0f, 0.0f, 1.0f}),
        RGBColor({0.5f, 0.0f, 1.0f}), RGBColor({1.0f, 0.0f, 1.0f}), RGBColor({1.0f, 0.0f, 0.5f}),
        RGBColor({0.6f, 0.0f, 0.0f}), RGBColor({0.6f, 0.3f, 0.0f}), RGBColor({0.6f, 0.6f, 0.0f}),
        RGBColor({0.3f, 0.6f, 0.0f}), RGBColor({0.0f, 0.6f, 0.0f}), RGBColor({0.0f, 0.6f, 0.3f}),
        RGBColor({0.0f, 0.6f, 0.6f}), RGBColor({0.0f, 0.3f, 0.6f}), RGBColor({0.0f, 0.0f, 0.6f}),
        RGBColor({0.3f, 0.0f, 0.6f}), RGBColor({0.6f, 0.0f, 0.6f}), RGBColor({0.6f, 0.0f, 0.3f})};
    constexpr auto numColors = sizeof(colors) / sizeof(colors[0]);

    return colors[index % numColors];
}
} // namespace PostscriptWriterColors

PostscriptWriter::PostscriptWriter(bool isTorus) : wrapAround(isTorus), ps_size{1020.0, 1020.0} {}

void PostscriptWriter::computeBoundaryBox(const std::vector<Point2D> &coordinates) {
    ps_min = {std::numeric_limits<coordinate>::max(), std::numeric_limits<coordinate>::max()};
    ps_max = {std::numeric_limits<coordinate>::min(), std::numeric_limits<coordinate>::min()};
    for (const auto &p : coordinates) {
        ps_min = ps_min.min(p);
        ps_max = ps_max.max(p);
    }

    ps_scale = (ps_size - (ps_border * 2.0)) / (ps_max - ps_min);
}

void PostscriptWriter::writeHeader(std::ofstream &file) const {
    /* Header */
    if (wrapAround) {
        file << "%!PS-Adobe-3.0 EPSF-3.0\n";
    } else {
        file << "%!PS-Adobe-1.0\n";
    }
    file << "%%Title: NetworKit visualization\n";
    file << "%%BoundingBox: 0.000 0.000 " << ps_size[0] << " " << ps_size[0] << "\n";
    file << "%%EndComments\n";
    if (!wrapAround) {
        file << "%%EndProlog\n";
        file << "gsave\n";
    }
}

void PostscriptWriter::writeMacros(std::ofstream &file) const {
    /* Macros */
    file << "/p {newpath} bind def\n"
            "/m {moveto} bind def\n"
            "/r {rmoveto} bind def\n"
            "/k {rlineto} bind def\n"
            "/l {lineto} bind def\n"
            "/n {rlineto} bind def\n"
            "/c {setrgbcolor} bind def\n"
            "/s {stroke} bind def\n"
            "/w {setlinewidth} bind def\n"
            "/h {show} bind def\n"
            "/a {arc closepath fill} bind def\n"
            "/b {closepath eofill} bind def\n";
}

void PostscriptWriter::writeClustering(const Graph &g, const std::vector<Point2D> &coordinates,
                                       const Partition &clustering, std::ofstream &file) {

    auto adjustToBoundingBox([&](Point2D p) { return (p - ps_min) * ps_scale + ps_border; });

    g.forEdges([&](node u, node v) {
        // set edge color
        if (clustering[u] == clustering[v] && clustering[u] != none) {
            // same cluster
            const auto color = PostscriptWriterColors::fromCyclicRotation(clustering[u]);
            file << color.r << " " << color.g << " " << color.b << " c ";
        } else {
            // different clusters -> grey
            file << "0.80 0.80 0.80 c 1.0 w ";
        }

        // set edge start and end point
        auto start = adjustToBoundingBox(coordinates[u]);
        auto end = adjustToBoundingBox(coordinates[v]);
        auto diff = end - start;

        if (wrapAround) {
            diff.apply([](index, coordinate val) { // TODO: externalize constants
                if (val > 500.0f)
                    return val - 1000.0f;
                if (val < -500.0f)
                    return val + 1000.0f;
                return val;
            });
            end = start + diff;
        }

        // write edge to file
        file << "p " << start[0] << " " << start[1] << " m " << end[0] << " " << end[1] << " l s\n";
    });

    // draw vertices
    float dotsize = 2.0;
    g.forNodes([&](node u) {
        if (clustering[u] != none) {
            // change color
            const auto color = PostscriptWriterColors::fromCyclicRotation(clustering[u]);
            file << color.r << " " << color.g << " " << color.b << " c ";
        } else {
            file << "0.0 0.0 0.0 c ";
        }

        const auto point = adjustToBoundingBox(coordinates[u]);

        file << "p " << point[0] << " " << point[1] << " " << dotsize << " 0.00 360.00 a s\n";
    });
}

void PostscriptWriter::init(std::ofstream &file) const {
    file.precision(3);
    file << std::fixed;

    writeHeader(file);
    writeMacros(file);
    file << "0.000 0.000 0.000 c\n";
}

void PostscriptWriter::write(const Graph &g, const std::vector<Point2D> &coordinates,
                             const Partition &clustering, const std::string &path) {
    assert(g.upperNodeIdBound() == coordinates.size());

    std::ofstream file(path);
    init(file);
    computeBoundaryBox(coordinates);

    writeClustering(g, coordinates, clustering, file);

    if (!wrapAround) {
        file << "grestore\n";
    }

    file.close();
}

void PostscriptWriter::write(const Graph &g, const std::vector<Point2D> &coordinates,
                             const std::string &path) {
    assert(g.upperNodeIdBound() == coordinates.size());

    ClusteringGenerator gen;
    Partition allNone = gen.makeOneClustering(g);
    write(g, coordinates, allNone, path);
}

} /* namespace NetworKit */
