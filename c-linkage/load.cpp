#include "cluster.hpp"


static bool is_whitespace(char c)
{
  if (c == ' ' or c == '\n' or c == '\t' or c == '\r') return true;
  return false;
}

static void split_string(vector<char*>* v, char* s, char sep_char)
{
  v->clear();
  char* t = s;
  bool in_word = false;
  while (*t != 0) {
    bool is_sep = false;
    if (sep_char > 0) {
      if (*t == sep_char)
        is_sep = true;
    }
    if (sep_char == 0) {
      if (is_whitespace(*t))
        is_sep = true;
    }
  
    if (not is_sep) {
      if (not in_word) v->push_back(t);
      in_word = true;
    } else {
      *t = 0;
      in_word = false;
    }
    t++;
  }
}

static inline int get_seq_index(char* seq)
{
  string s(seq);
  if (map_name_to_index.find(s) == map_name_to_index.end()){
    fprintf(stderr, "error. '%s' not found in names file.\n", seq);
    exit(1);
  }
  return map_name_to_index[s];
}

#define PERM_BUF_SIZE 1024*1024
char perm_buf[PERM_BUF_SIZE];
char* perm_buf_start = perm_buf;
char* perm_buf_end = perm_buf;
static char* read_line(char* buffer, size_t buf_size, FILE* f){
  size_t l = perm_buf_end - perm_buf_start;
  if (l < buf_size) {
    if (perm_buf_end - perm_buf_start > 0){
      memmove(perm_buf, perm_buf_start, l);
    }
    perm_buf_start = perm_buf;
    perm_buf_end = perm_buf + l;
    size_t r = fread(perm_buf_end, 1, PERM_BUF_SIZE - l, f);
    perm_buf_end += r;
    //fprintf(stderr, "  read %d\n", r);
  }
  //fprintf(stderr, "  %p,  %p\n", perm_buf_start, perm_buf_end);
  if (perm_buf_start >= perm_buf_end) return NULL;
  
  char* c = perm_buf_start;
  while (*c != '\n' and c != perm_buf_end) c++;
  *c = 0;
  char* t = perm_buf_start;
  perm_buf_start = c + 1;
  strcpy(buffer, t);
  return t;
}

#define LINE_SZ 100000

void load_names(const char* names_fn)
{
  FILE* f = fopen(names_fn, "r");

  char buffer[LINE_SZ];
  el_index seq_index = 0;
  vector<char*> words;
  
  while(fgets(buffer, LINE_SZ, f)){
    split_string(&words, buffer, 0);
    if (words.size() == 0) continue;
    
    if (words.size() != 2) {
      fprintf(stderr, "error no 2 columns in %s\n", names_fn);
      fprintf(stderr, "output of columns follows now\n");
      for (int i = 0; i < int(words.size()); i++){
        fprintf(stderr, "%d: %s\n", i, words[i]);
      }
      exit(1);
    }
    
    // Split the second word across commas
    vector<char*> seqs;
    split_string(&seqs, words[1], ',');

    // Create vector of strings
    vector<string> list;
    for (int i = 0; i < int(seqs.size()); i++){
      string s(seqs[i]);
      list.push_back(s);
    }
    assert(list.size() > 0);
    map_names[seq_index] = list;

    // Form a string of first word, and map it.
    string s(words[0]);
    map_name_to_index[s] = seq_index;
    map_index_to_mult[seq_index] = seqs.size();
    assert(seqs.size() > 0);
    num_total += seqs.size();
    seq_index++;
    
  }
  num_uniq = seq_index;
  num_el = num_uniq * (num_uniq - 1) / 2;
  fclose(f);
}

void load_dist(const char* dist_fn)
{
  FILE* f = fopen(dist_fn, "r");
  
  char buffer[LINE_SZ];
  dist_index r_dists = 0;
  vector<char*> words;
  
  while (read_line(buffer, 1000, f)){
    split_string(&words, buffer, 0);
    // empty line.
    if (words.size() == 0) continue;
    
    if (words.size() != 3) {
      fprintf(stderr, "error. not 3 columns in %s\n", dist_fn);
      fprintf(stderr, "output of columns follows now\n");
      for (int i = 0; i < int(words.size()); i++){
        printf("%d: %s\n", i, words[i]);
      }
      exit(1);    
    }
    
    el_index i1 = get_seq_index(words[0]);
    el_index i2 = get_seq_index(words[1]);
    dist d = atof(words[2]);
    assert(i1 != i2);
    set_dist(i1, i2, d);

    r_dists++;
    dist_index step_size = num_el / 10;
    if (r_dists % step_size == 0)
      fprintf(stderr, "  read ~%" PRI_dist_index "%%\n", r_dists * 100 / num_el);
  }
  fclose(f);  
  fprintf(stderr, "  read %" PRI_dist_index " distances and expected %" PRI_dist_index "!\n", 
          r_dists, num_el);
  if (r_dists != num_el){
    exit(1);
  }
}

