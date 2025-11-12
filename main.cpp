
#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
#include <openvdb/openvdb.h>
#include <tira/volume.h>
#include <openvdb/tools/ChangeBackground.h>

#include <openvdb/tools/Morphology.h>

// input arguments
std::string in_inputname;
std::string in_outputname;
std::string in_templates;
int in_cuda;
unsigned int in_foreground;

struct coord {
    int dx;
    int dy;
    int dz;
    int value;
};

struct Template {
    char label;
    unsigned int id;
    std::vector<coord> coords;
};

struct template_entry {
    char label;
    unsigned int id;
    int x;
    int y;
    int z;
    int value;
};

template_entry process_entry(std::string line, int line_num) {
    template_entry e;
    std::stringstream ss(line);                     // create a string stream object to parse the line
    char comma;
    if (!(ss >> e.label >> comma >> e.id >> comma >> e.x >> comma >> e.y >> comma >> e.z >> comma >> e.value)) {
        throw std::runtime_error("Could not parse line " + line_num);
    }
    return e;
}

// load templates from CSV
std::vector<Template> load_templates_from_csv(const std::string& filename) {
    std::ifstream file(filename);                                                   // open the template csv file for reading
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }

    std::vector<Template> templates;                                                // create an empty list of templates
    //std::map<std::pair<char, int>, Template> template_map; // to store (class, ID) -> Template



    std::string line;                               // create a line used to buffer a line of text in the input file
    std::getline(file, line);                   // throw away the first line

    unsigned int line_number = 0;
    unsigned int template_number = 0;
    Template t;                                         // temporary template used to store values while reading
    while (std::getline(file, line)) {              // while there are still lines in the file, read a line
        line_number++;                                  // increment the line number

        template_entry e = process_entry(line, line_number);    // process an entry from the CSV file
        if (line_number == 1) {                                     // if this is the first line
            t.label = e.label;                                      // create the first template
            t.id = e.id;
            t.coords.clear();
        }
        else if (t.label != e.label || t.id != e.id) {              // if the current line represents a new template
            templates.push_back(t);                                 // store the previous template
            t.label = e.label;                                      // update the label and ID
            t.id = e.id;
            t.coords.clear();                                       // empty the coordinates
        }
        coord c;                                                    // create and fill the coordinate represented by this line
        c.dx = e.x;
        c.dy = e.y;
        c.dz = e.z;
        c.value = e.value;  //  assign the value field!
        t.coords.push_back(c);                                      // add the coordinate to the template

    }

    // map values to vector
    //for (const auto& pair : template_map) {
    //    templates.push_back(pair.second);
    //}

    return templates;
}



struct Point {
    int x, y, z;
};

// Initialize a zero OpenVDB grid 
openvdb::FloatGrid::Ptr zero_vdb_volume()
{
    // Background is 0.0f and no voxels active by default
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);
    grid->setTransform(openvdb::math::Transform::createLinearTransform(1.0)); // optional voxel size
    return grid;
}



//  is_tail_point (endpoint detector)

bool is_tail_point(openvdb::FloatGrid::Accessor& acc, int x, int y, int z)
{
    static const Point n26[] = {
        { 1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
        { 1,1,0},{ 1,-1,0},{-1,1,0},{-1,-1,0},
        { 1,0,1},{ 1,0,-1},{-1,0,1},{-1,0,-1},
        { 0,1,1},{ 0,1,-1},{ 0,-1,1},{ 0,-1,-1},
        { 1,1,1},{ 1,1,-1},{ 1,-1,1},{ 1,-1,-1},
        {-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}
    };

    int neighbor_count = 0;

    for (const auto& o : n26) {
        openvdb::Coord c(x + o.x, y + o.y, z + o.z);
        const float v = acc.getValue(c); // inactive voxels return background
        if (v > 127.5f) ++neighbor_count;
    }

    // endpoint or near-endpoint
    return (neighbor_count == 1 || neighbor_count == 2);
}



// matches_template (template pattern check)

bool matches_template(openvdb::FloatGrid::Accessor& acc,
    int x, int y, int z, const Template& tmpl)
{
    bool has_at_least_one_1 = false;
    int x_count = 0;

    for (const auto& coord : tmpl.coords) {
        const int nx = x + coord.dx;
        const int ny = y + coord.dy;
        const int nz = z + coord.dz;
        const openvdb::Coord c(nx, ny, nz);
        const float v = acc.getValue(c);

        if (coord.value == -1) { // wildcard X
            ++x_count;
            if (v > 127.5f) has_at_least_one_1 = true;
        }
        else if (std::abs(v - float(coord.value)) > 1e-3f) {
            return false;
        }
    }

    return (x_count > 0) ? has_at_least_one_1 : true;
}



// is_p_simple 

bool is_p_simple(openvdb::FloatGrid::Accessor& acc, int x, int y, int z)
{
    int T26 = 0, T6 = 0;
    std::vector<Point> R_neighbors;
    R_neighbors.reserve(26);

    // Count neighbors in 26-neighborhood
    for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
            {
                if (!dx && !dy && !dz) continue;
                const openvdb::Coord c(x + dx, y + dy, z + dz);
                const float v = acc.getValue(c);

                if (v > 127.5f) {
                    ++T26;
                }
                else {
                    R_neighbors.push_back({ x + dx, y + dy, z + dz });
                }
            }

    // 6-neighborhood background count
    static const Point n6[] = {
        { 1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}
    };
    for (const auto& o : n6) {
        const openvdb::Coord c(x + o.x, y + o.y, z + o.z);
        const float v = acc.getValue(c);
        if (v < 127.5f) ++T6;
    }

    // Adjacency check among background neighbors
    for (const auto& a : R_neighbors) {
        bool found = false;
        for (const auto& b : R_neighbors) {
            if (std::abs(b.x - a.x) <= 1 &&
                std::abs(b.y - a.y) <= 1 &&
                std::abs(b.z - a.z) <= 1) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }

    // P-simple condition
    return (T26 > 1 && T6 > 1);
}


#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>

openvdb::FloatGrid::Ptr thinning_lohou_parallel(
    openvdb::FloatGrid::Ptr& volume,
    const std::vector<Template>& templates,
    int num_threads = tbb::info::default_concurrency())
{
    using Grid = openvdb::FloatGrid;
    using Tree = Grid::TreeType;

    openvdb::FloatGrid::Ptr output = volume->deepCopy();
    output->setTransform(volume->transform().copy());

    bool changed = true;
    int iter = 0;

    // Limit parallelism
    tbb::global_control gc(tbb::global_control::max_allowed_parallelism, num_threads);

    while (changed && iter < 200) {
        changed = false;
        std::atomic<size_t> numDeleted(0);

        openvdb::FloatGrid::Ptr to_delete = openvdb::FloatGrid::create(0.0f);
        to_delete->setTransform(output->transform().copy());

        using LeafNode = Tree::LeafNodeType;
        using LeafIterRange = openvdb::tree::IteratorRange<Tree::LeafCIter>;

        //Parallel deletion marking
        tbb::parallel_for(
            LeafIterRange(output->tree().cbeginLeaf()),
            [&](LeafIterRange& range) {
                auto acc = output->getAccessor();   // thread-local accessor
                auto delAcc = to_delete->getAccessor(); //  thread-local accessor

                for (; range; ++range) {
                    auto& leaf = const_cast<LeafNode&>(*range.iterator());
                    for (auto iter = leaf.beginValueOn(); iter; ++iter) {
                        const openvdb::Coord c = iter.getCoord();
                        const float v = iter.getValue();
                        if (std::abs(v - 255.0f) > 1e-3f) continue;

                        const int x = c.x(), y = c.y(), z = c.z();

                        if (!is_tail_point(acc, x, y, z)) {
                            for (const auto& tmpl : templates) {
                                if (matches_template(acc, x, y, z, tmpl) &&
                                    is_p_simple(acc, x, y, z)) {
                                    delAcc.setValueOn(c, 255.0f);
                                    numDeleted.fetch_add(1, std::memory_order_relaxed);
                                    break;
                                }
                            }
                        }
                    }
                }
            });

        // deletions in parallel 
        if (numDeleted > 0) {
            changed = true;
            tbb::parallel_for(
                LeafIterRange(to_delete->tree().cbeginLeaf()),
                [&](LeafIterRange& range) {
                    auto acc = output->getAccessor(); // thread-local accessor
                    for (; range; ++range) {
                        auto& leaf = const_cast<LeafNode&>(*range.iterator());
                        for (auto iter = leaf.beginValueOn(); iter; ++iter) {
                            acc.setValue(iter.getCoord(), 0.0f);
                        }
                    }
                });
            output->tree().prune();
        }

        ++iter;
        std::cout << "[iter " << iter << "] deleted " << numDeleted
            << " active: " << output->activeVoxelCount() << std::endl;
    }

    std::cout << "Total iterations: " << iter << std::endl;
    return output;
}





// Dilate active voxels by N layers and then fill band with 0 values
openvdb::FloatGrid::Ptr add_zero_band(openvdb::FloatGrid::Ptr grid, int radius = 3)
{
    using namespace openvdb::tools;


    openvdb::FloatGrid::Ptr dilated = openvdb::FloatGrid::create(0.0f);
    dilated->setTransform(grid->transform().copy());
    dilated->tree().topologyUnion(grid->tree());

    dilateActiveValues(dilated->tree(), radius);


    auto inAcc = grid->getAccessor();
    auto outAcc = dilated->getAccessor();

    for (auto iter = dilated->beginValueOn(); iter; ++iter) {
        openvdb::Coord c = iter.getCoord();
        if (inAcc.isValueOn(c)) {
            outAcc.setValue(c, 255.0f); // foreground
        }
        else {
            outAcc.setValue(c, 0.0f);   // background 
        }
    }

    dilated->tree().prune(); // clean up empty tiles
    return dilated;
}


//  already-active voxels, background stays inactive.
template <typename GridPtr>
void binarize_active_only(GridPtr grid, float threshold = 0.5f)
{
    using Grid = typename GridPtr::element_type;
    typename Grid::Accessor acc = grid->getAccessor();

    for (auto it = grid->beginValueOn(); it; ++it) {
        const float v = it.getValue();
        if (v >= threshold) {
            it.setValue(255.0f);   // Foreground
        }
        else {
            it.setValueOff();      // Background - deactivate voxel
        }
    }

    grid->tree().prune();          // Remove empty tiles
}


int main() {
    openvdb::initialize();

    std::vector<Template> templates =
        load_templates_from_csv("templates_255.csv");

    std::string input_filename = "cube_500.vdb";
    openvdb::io::File file_sdf(input_filename);
    file_sdf.open();
    openvdb::GridBase::Ptr baseGrid =
        file_sdf.readGrid(file_sdf.beginName().gridName());
    file_sdf.close();

    openvdb::FloatGrid::Ptr input_grid =
        openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

    binarize_active_only(input_grid, 0.5f);
    input_grid = add_zero_band(input_grid, 3);

    openvdb::tools::changeBackground(input_grid->tree(), 0);
    input_grid->tree().prune();
    input_grid->pruneGrid(0.0f);
    openvdb::tools::dilateActiveValues(input_grid->tree(), 1);

    // different thread counts
    for (int threads : {1,2,4,8,16}) {
        auto start = std::chrono::high_resolution_clock::now();
        auto thinned = thinning_lohou_parallel(input_grid, templates, threads);
        auto end = std::chrono::high_resolution_clock::now();

        double elapsed = std::chrono::duration<double>(end - start).count();
        std::cout << "[Threads " << threads << "] Time: " << elapsed << " s\n";

        
        std::string output_filename = "thin_lohou_500.vdb";
        openvdb::io::File output_file(output_filename);
        openvdb::GridPtrVec grids;
        grids.push_back(thinned);

        try {
            output_file.write(grids);
            output_file.close();
            std::cout << "Successfully saved thinned VDB to: " << output_filename << std::endl;
        }
        catch (const std::exception& e) {
            std::cerr << "Error saving VDB file: " << e.what() << std::endl;
            return 1;
        }

       
    }

    return 0;
}


