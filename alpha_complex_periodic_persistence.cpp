#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <random>
#include <cmath>

#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>


using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = Simplex_tree::Filtration_value;
using Persistent_cohomology =
    Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Gudhi::persistent_cohomology::Field_Zp>;


template <typename AlphaComplex3d>
std::vector<typename AlphaComplex3d::Point_3> coords_to_points(const std::vector<std::vector<double>> &coords) {
    using Point = typename AlphaComplex3d::Point_3;

    std::vector<Point> points;
    for (int i = 0; i < coords.size(); i++) {
        points.push_back(Point(coords[i][0], coords[i][1], coords[i][2]));
    }

    return points;
}


std::vector<std::vector<std::vector<double>>> calc_persistence(const std::vector<std::vector<double>> &coords, int coeff_field, double min_persistence, double boxsize = 0.) {

    Simplex_tree simplex_tree;

    if (boxsize > 0) {
        using Alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, true>;
        auto points = coords_to_points<Alpha_complex_3d>(coords);
        Alpha_complex_3d alpha_complex(points, -boxsize, -boxsize, -boxsize, boxsize, boxsize, boxsize);
        alpha_complex.create_complex(simplex_tree);
    } else {
        using Alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, false>;
        auto points = coords_to_points<Alpha_complex_3d>(coords);
        Alpha_complex_3d alpha_complex(points);
        alpha_complex.create_complex(simplex_tree);
    }

    int st_dimension = simplex_tree.dimension();

    Persistent_cohomology pcoh(simplex_tree, true);
    pcoh.init_coefficients(coeff_field);
    pcoh.compute_persistent_cohomology(min_persistence);

    std::vector<std::vector<std::vector<double>>> persistence_pairs;
    std::vector<std::pair<Filtration_value, Filtration_value>> pairs_d;
    std::vector<std::vector<double>> pairs_d_vec;
    for (int dim = 0; dim < st_dimension; dim++) {
        pairs_d = pcoh.intervals_in_dimension(dim);
        pairs_d_vec.clear();
        for (int i = 0; i < pairs_d.size(); i++) {
            pairs_d_vec.push_back({sqrt(pairs_d[i].first), sqrt(pairs_d[i].second)});
        }
        persistence_pairs.push_back(pairs_d_vec);
    }

    return persistence_pairs;
}


int main() {
    int N_points = 10;
    double boxsize = 1;

    std::default_random_engine generator;
    generator.seed(10);
    std::uniform_real_distribution<double> distribution(-boxsize, boxsize);

    std::vector<std::vector<double>> points;
    for (int i = 0; i < N_points; i++) {
        std::vector<double> point = {distribution(generator), distribution(generator), distribution(generator)};
        points.push_back(point);
    }

    std::cout << "Points:" << std::endl;
    for (int i = 0; i < N_points; i++) {
        std::cout << "\t" << points[i][0] << ", " << points[i][1] << ", " << points[i][2] << std::endl;
    }

    std::vector<std::vector<std::vector<double>>> persistence_pairs = calc_persistence(points, 2, 0., boxsize);

        for (int dim = 0; dim < persistence_pairs.size(); dim++) {
        std::cout << "H" << dim << std::endl;
        for (int i = 0; i < persistence_pairs[dim].size(); i++) {
            std::cout << "\t" << persistence_pairs[dim][i][0] << "  " << persistence_pairs[dim][i][1] << std::endl;
        }
    }
}
