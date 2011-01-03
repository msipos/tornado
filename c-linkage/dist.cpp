#include "cluster.hpp"

dist* dist_matrix;
el_index* sort_r;
el_index* sort_c;

// assuming r < c.
#define DIST(r, c)  (dist_index(r) * dist_index(num_uniq) + dist_index(c) - (dist_index(r) + 1)*(dist_index(r) + 2)/2)

void init_dist()
{
  fprintf(stderr, "  using ~%" PRI_dist_index " Mb of RAM\n", (sizeof(float) + 2*sizeof(int))
                                                              * num_el / 1024 / 1024);
  // Initialize distance matrix
  dist_matrix = (dist*) malloc(sizeof(dist) * num_el);

  // Allocate sorted element pointers
  sort_r = (el_index*) malloc(sizeof(el_index) * num_el);
  sort_c = (el_index*) malloc(sizeof(el_index) * num_el);
  
  // Initialize sort_r and sort_c
  for (el_index r = 0; r < num_uniq; r++){
    for (el_index c = r + 1; c < num_uniq; c++){
      sort_r[DIST(r,c)] = r;
      sort_c[DIST(r,c)] = c;
    }
  }
}

void set_dist(el_index a, el_index b, dist d)
{
  if (a == b) return;
  if (a < b) {
    dist_matrix[DIST(a,b)] = d;
  } else {
    dist_matrix[DIST(b,a)] = d;
  }
}

dist get_dist(el_index a, el_index b)
{
  if (a == b) return 0.0;
  if (a < b) {
    return dist_matrix[DIST(a,b)];
  } else {
    return dist_matrix[DIST(b,a)];
  }
}

void display_matrix()
{
  for (el_index r = 0; r < num_uniq; r++){
    for (el_index c = 0; c < num_uniq; c++){
      float f = get_dist(r,c);
      if (f >= 0.0) {
        fprintf(stderr, " %.3f ", get_dist(r, c));
      } else {
        fprintf(stderr, "%.3f ", get_dist(r, c));
      }
    }
    fprintf(stderr, "\n");
  }
}

