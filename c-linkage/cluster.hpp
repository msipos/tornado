#include <vector>
#include <utility>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

using namespace std;

// Index into the distance matrix.
typedef uint64_t dist_index;
#define PRI_dist_index PRIu64

// Index of an element.
typedef uint32_t el_index;
#define PRI_el_index PRIu32

// Type of a distance.
typedef float dist;
#define PRI_dist "f"

// Maps name of an element to its index.
extern map <string, el_index> map_name_to_index;
// Maps index of an element to its multiplicity (i.e. how many elements stand for this element).
extern map <el_index, el_index> map_index_to_mult;
// Names file map: maps index of a reference sequence to vector of sequence names
extern map <el_index, vector<string> > map_names;
// Basically map above, but faster.
extern el_index* vec_index_to_mult;
// Number of unique elements.
extern el_index num_uniq;
// Total number of elements, i.e. (# of unique elements) * (average multiplicity).
extern el_index num_total;
// Number of matrix elements.
extern dist_index num_el;
// Distance matrix.
extern dist* dist_matrix;
// Sorted distance matrix r,c pairs
extern el_index* sort_r;
extern el_index* sort_c;

// load.cpp
void load_names(const char* names_fn);
void load_dist(const char* dist_fn);

// dist.cpp
void set_dist(el_index a, el_index b, dist d);
dist get_dist(el_index a, el_index b);
void display_matrix();
void init_dist();
