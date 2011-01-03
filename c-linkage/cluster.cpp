#include "cluster.hpp"

map <string, el_index> map_name_to_index;
map <el_index, el_index> map_index_to_mult;
map <el_index, vector<string> > map_names;
el_index* vec_index_to_mult;
el_index num_uniq;
el_index num_total;
dist_index num_el;

// Clusters array
vector<el_index>** clusters;

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

static inline void qsort_swap(dist_index a, dist_index b)
{
  if (a == b) return;

  el_index t_r, t_c;
  t_r = sort_r[a];
  t_c = sort_c[a];
  sort_r[a] = sort_r[b];
  sort_c[a] = sort_c[b];
  sort_r[b] = t_r;
  sort_c[b] = t_c;
}

#define QSDIST(a)  get_dist(sort_r[a], sort_c[a])

static inline bool qsort_gt(dist_index a, dist_index b)
{
  return QSDIST(a) > QSDIST(b);
}

dist_index qsorted = 0;
dist_index prev_perc = 0;
static inline void qsort_update()
{
  dist_index new_perc = qsorted * 100 / num_el;
  if (new_perc > prev_perc and new_perc % 5 == 0){
    fprintf(stderr,"  sorted ~%" PRI_dist_index "%%\n", new_perc);
  }
  prev_perc = new_perc;
}

// Sort [x1,x2] portion of the sort_r and sort_c arrays
void qsort_r(dist_index x1, dist_index x2)
{
  //fprintf(stderr, "  qsort_r(%d, %d)\n", x1, x2);

  dist_index l = x2 - x1 + 1;
  assert(l >= 0);
  if (l <= 1) {
    qsorted += l;
    qsort_update();
    return;
  }
  if (l == 2) {
    qsorted += 2;
    qsort_update();
    if (qsort_gt(x1, x2)){
      qsort_swap(x1, x2);
    }
    return;
  }

  // Find the pivot.
  dist_index xm = (x1 + x2) / 2;
  dist d1 = QSDIST(x1);
  dist d2 = QSDIST(x2);
  dist dm = QSDIST(xm);
  dist_index xp;
  if (dm > d1) {
    if (d2 > dm) {
      // d2 > dm > d1
      xp = xm;
    } else {
      if (d2 > d1) {
        // dm > d2 > d1
        xp = x2;
      } else {
        // dm > d1 > d2
        xp = x1;
      }
    }
  } else {
    // d1 > dm
    if (dm > d2) {
      // d1 > dm > d2
      xp = xm;
    } else {
      // d2 > dm
      if (d1 > d2) {
        // d1 > d2 > dm
        xp = x2;
      } else {
        // d2 > d1 > dm
        xp = x1;
      }
    }
  }

  // Make x2 be the pivot.
  qsort_swap(x2, xp);

  dist_index t = x1;
  for (dist_index n = x1; n < x2; n++){
    if (not qsort_gt(n, x2)) {
      qsort_swap(t, n);
      t++;
    }
  }
  qsort_swap(t, x2);
  qsorted += 1;
  qsort_update();

  qsort_r(x1, t - 1);
  qsort_r(t + 1, x2);
}

// Temporary storage.
el_index* tmp_r;
el_index* tmp_c;

static inline void merge_set(dist_index dest, dist_index src)
{
  tmp_r[dest] = sort_r[src];
  tmp_c[dest] = sort_c[src];
}

void merge_sort(dist_index x1, dist_index x2)
{
  dist_index l = x2 - x1 + 1;
  if (l <= 1) {
    qsorted += l;
    qsort_update();
    return;
  }
  dist_index xm = (x1 + x2) / 2;
  merge_sort(x1, xm);
  merge_sort(xm+1, x2);

  // Finally merging step.
  dist_index n1 = x1;
  dist_index n2 = xm+1;

  for(dist_index n = x1; n <= x2; n++){
//    fprintf(stderr,
//            "n1 = %" PRI_dist_index ", n2 = %" PRI_dist_index ", n = %" PRI_dist_index "\n"
//            "d1 = %" PRI_dist ", d2 = %" PRI_dist "\n", n1, n2, n, QSDIST(n1), QSDIST(n2));
    // n1 blew through its section
    if (n1 > xm){
      assert(n2 <= n);
      merge_set(n, n2);
      n2++;
      continue;
    }

    // n2 blew through its section
    if (n2 > x2){
      assert(n1 <= n);
      merge_set(n, n1);
      n1++;
      continue;
    }

    // or both n1 and n2 are valid
    if (QSDIST(n1) > QSDIST(n2)){
      merge_set(n, n2);
      n2++;
    } else {
      merge_set(n, n1);
      n1++;
    }
  }

  for(dist_index n = x1; n <= x2; n++){
    sort_r[n] = tmp_r[n];
    sort_c[n] = tmp_c[n];
  }

//  fprintf(stderr, "  \n");
//  dist f = -1.0;
//  for (el_index i = x1; i <= x2; i++){
//    assert(QSDIST(i) >= f);
//    fprintf(stderr, "  %f\n", QSDIST(i));
//    f = QSDIST(i);
//  }
//  fprintf(stderr, "  \n");
}

void init_sort()
{
  // Allocate temporary storage
  tmp_r = (el_index*) malloc(num_el * sizeof(el_index));
  tmp_c = (el_index*) malloc(num_el * sizeof(el_index));

  merge_sort(0, num_el - 1);
//  for (el_index i = 0; i < num_el; i++){
//    fprintf(stderr, "  %f\n", QSDIST(i));
//  }
  free(tmp_r);
  free(tmp_c);

  // Sort distance matrix elements
  //qsort_r(0, num_el - 1);
}

#define MARKER (-1.0)

void merge(el_index a, el_index b)
{
  vector<el_index>* c_a = clusters[a];
  vector<el_index>* c_b = clusters[b];
  el_index sz_a = c_a->size();
  el_index sz_b = c_b->size();

  // Destroy all links between clusters[a] and clusters[b]
  for (el_index i = 0; i < sz_a; i++){
    for (el_index j = 0; j < sz_b; j++){
      set_dist(c_a->at(i), c_b->at(j), MARKER);
    }
  }

  // Preserve only maximum links to e
  for (el_index e = 0; e < num_uniq; e++){
    dist d1 = MARKER;
    el_index e1 = 0;
    dist d2 = MARKER;
    el_index e2 = 0;

    for (el_index i = 0; i < sz_a; i++){
      e1 = c_a->at(i);
      if (e1 == e) goto end;
      d1 = get_dist(e1, e);
      if (d1 >= 0.0) break;
    }

    for (el_index i = 0; i < sz_b; i++){
      e2 = c_b->at(i);
      if (e2 == e) goto end;
      d2 = get_dist(e2, e);
      if (d2 >= 0.0) break;
    }

    if (d1 > d2)
      set_dist(e2, e, MARKER);
    else
      set_dist(e1, e, MARKER);

    end:
    ;
  }

  // Now merge the two clusters.
  for (el_index j = 0; j < sz_b; j++){
    el_index e = c_b->at(j);
    c_a->push_back(e);
    clusters[e] = c_a;
  }

  // Now, go through each cluster that exists, and make sure only one
  // distance exists to it.
  set<vector<el_index>*> s;
  for (el_index i = 0; i < num_uniq; i++){
    if (s.find(clusters[i]) == s.end()){
      s.insert(clusters[i]);
      vector<el_index>* cl = clusters[i];
      el_index sz = cl->size();

      // Find maximum distance
      dist dmax = MARKER;
      for (el_index j = 0; j < sz; j++){
        el_index e1 = cl->at(j);
        for (el_index k = 0; k < sz_a + sz_b; k++){
          el_index e2 = c_a->at(k);
          dist d = get_dist(e1, e2);
          if (d > dmax) dmax = d;
        }
      }

      //fprintf(stderr, "  found dmax = %f\n", dmax);

      for (el_index j = 0; j < sz; j++){
        el_index e1 = cl->at(j);
        for (el_index k = 0; k < sz_a + sz_b; k++){
          el_index e2 = c_a->at(k);
          dist d = get_dist(e1, e2);
          if (d < dmax and d >= 0.0) {
            set_dist(e1, e2, MARKER);
          }
        }
      }
    }
  }
}

void do_rabund(float R)
{
  vector<el_index> rabund;
  vector<string> list;
  set<vector<el_index>*> s;
  string single;

  for (el_index i = 0; i < num_uniq; i++){
    if (s.find(clusters[i]) == s.end()){
      s.insert(clusters[i]);
      single = "";
      vector<string> temp;
      for (int j = 0; j < int(clusters[i]->size()); j++){
        el_index e = clusters[i]->at(j);
        vector<string> temp2 = map_names[e];
        assert(temp2.size() > 0);
        for (int k = 0; k < int(temp2.size()); k++){
          temp.push_back(temp2[k]);
        }
      }
      assert(temp.size() > 0);
      single += temp[0];
      for (int j = 1; j < int(temp.size()); j++){
        single += ",";
        single += temp[j];
      }
      list.push_back(single);

      // For rabund:
      int total = 0;
      for (int j = 0; j < int(clusters[i]->size()); j++){
        total += vec_index_to_mult[clusters[i]->at(j)];
      }
      rabund.push_back(total);
    }
  }
  sort(rabund.begin(), rabund.end());
  reverse(rabund.begin(), rabund.end());

  FILE* fl = fopen("cluster.list", "a");
  if (R < 0.0)
    fprintf(fl, "unique %" PRI_el_index " ", el_index(rabund.size()));
  else
    fprintf(fl, "%.6f   %" PRI_el_index " ", R, el_index(rabund.size()));

  for (int i = 0; i < int(list.size()); i++){
    fprintf(fl, "%s ", list[i].c_str());
  }
  fprintf(fl, "\n");
  fclose(fl);

  fprintf(stderr, "  OTUs left to merge: %u\n", rabund.size());

  FILE* f = fopen("cluster.rabund", "a");
  if (R < 0.0)
    fprintf(f, "unique %" PRI_el_index " ", el_index(rabund.size()));
  else
    fprintf(f, "%.6f   %" PRI_el_index " ", R, el_index(rabund.size()));
  for (int i = 0; i < int(rabund.size()); i++){
    fprintf(f, "%d ", rabund[i]);
  }
  fprintf(f, "\n");
  fclose(f);
}

#define MODE_EVERY  1
#define MODE_ONLY   2

int mode;
double maxdist;
int every;
double only;

void cluster()
{
  if (mode == MODE_EVERY){
    do_rabund(-1.0);
  }
  el_index i = 0;

  for (el_index n = 0; n < num_el; n++){
    dist R = QSDIST(n);
    if (R < 0.0) continue;

    if (R > maxdist) break;

    el_index r = sort_r[n];
    el_index c = sort_c[n];

    merge(r, c);
    if (mode == MODE_EVERY){
      if (i % every == 0) do_rabund(R);
      i++;
    }

    if (mode == MODE_ONLY){
      double difference = R - only;
      if (difference < 0.0) difference *= -1.0;
      if (R > only or difference / only < 0.01){
        do_rabund(R);
        return;
      }
    }
  }
}

void verify_sortedness()
{
  dist f = 0.0;
  for (el_index i = 0; i < num_el; i++){
    if (QSDIST(i) < f)
      fprintf(stderr, "  Error with sorting!\n");
    f = QSDIST(i);
  }
}

void init()
{
  // Initialize vec_index_to_mult
  vec_index_to_mult = (el_index*) malloc(sizeof(el_index) * num_uniq);
  for (el_index i = 0; i < num_uniq; i++){
    vec_index_to_mult[i] = map_index_to_mult[i];
  }

  // Allocate clusters array
  clusters = (vector<el_index>**) malloc(sizeof(vector<el_index>*) * num_uniq);
  for (el_index i = 0; i < num_uniq; i++){
    clusters[i] = new vector<el_index>;
    clusters[i]->push_back(i);
  }
}

void usage()
{
  fprintf(stdout,
          "Usage:  \n"
          "  c-linkage DIST_FILE NAMES_FILE [MODIFIERS...]\n"
          "\n"
          "where DIST_FILE is a distance matrix in the \"column\" format. E.g.:\n"
          "    SEQ_ID_1 SEQ_ID_2 0.3\n"
          "    SEQ_ID_2 SEQ_ID_3 0.4\n"
          "    SEQ_ID_1 SEQ_ID 3 0.2\n"
          "and NAMES_FILE is a file that describes the unique sequences in DIST_FILE.\n"
          "That is, NAMES_FILE indicates which sequences are identical to those in\n"
          "DIST_FILE. For example,\n"
          "    SEQ_ID_1 SEQ_ID_1\n"
          "    SEQ_ID_2 SEQ_ID_2\n"
          "    SEQ_ID_3 SEQ_ID_3,SEQ_ID_4\n"
          "In the example above, SEQ_ID_1 and SEQ_ID_2 are unique sequences, whereas\n"
          "SEQ_ID_3 and SEQ_ID_4 are identical to each other.  SEQ_ID_3 and SEQ_ID_4\n"
          "are represented under SEQ_ID_3's name in the distance matrix.(Note: second\n"
          "column of NAMES_FILE must be comma delimited.)\n");
}

int main(int argc, char** argv)
{
  if (argc < 3) {
    usage();
    return 0;
  }
  every = 10;
  maxdist = 1.0;
  mode = MODE_EVERY;

  for (int i = 3; i < argc; i+=2){
    char* arg = argv[i];

    if (i+1 == argc){
      usage();
      return 0;
    }

    if (strcmp(arg, "upto") == 0){
      maxdist = atof(argv[i+1]);
    } else if (strcmp(arg, "every") == 0){
      every = atoi(argv[i+1]);
      mode = MODE_EVERY;
    } else if (strcmp(arg, "only") == 0){
      only = atof(argv[i+1]);
      mode = MODE_ONLY;
    } else {
      usage();
      return 0;
    }
  }

  char* dist_fn = argv[1];
  char* names_fn = argv[2];

  fprintf(stderr, "loading names file '%s'...\n", names_fn);
  load_names(names_fn);
  fprintf(stderr, "according to names file:\n  number of unique seqs: %" PRI_el_index "\n"
                  "  total seqs: %" PRI_el_index "\n", num_uniq, num_total);
  fprintf(stderr, "allocating matrix...\n");
  init_dist();
  init();
  fprintf(stderr, "loading distance matrix file '%s'...\n", dist_fn);
  load_dist(dist_fn);
  fprintf(stderr, "sorting distance matrix...\n");
  init_sort();
  verify_sortedness();
  fprintf(stderr, "clustering...\n");

  // Erase the two output files.
  FILE* f = fopen("cluster.rabund", "w");
  fclose(f);
  f = fopen("cluster.list", "w");
  fclose(f);

  cluster();
}
