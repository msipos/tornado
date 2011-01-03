#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cassert>

using namespace std;

vector<char*> seqs;
vector<char*> names;

static void usage(char* p)
{
  fprintf(stderr, "\n"
                  "Usage:\n"
                  "   %s FASTA_FILE  [ DIST_MATRIX_FILE ]\n"
                  "\n"
                  "  If DIST_MATRIX_FILE is not given, this program will\n"
                  "  calculate the PSI distance matrix of the sequences given\n"
                  "  in FASTA_FILE.\n"
                  "\n"
                  "  If DIST_MATRIX_FILE is given, this program will output\n"
                  "  the distance matrix from DIST_MATRIX_FILE *normalized*\n"
                  "  in a way such that the average distance in the output\n"
                  "  will be equal to the average distance in the PSI distance\n"
                  "  matrix.\n"
                  "\n", p);
}

static void add_seq(char* name, char* seq)
{
  seqs.push_back(strdup(seq));
  names.push_back(strdup(name));
}

static void load_fasta(char* filename)
{
  // Load FASTA file
  FILE* f = fopen(filename, "r");
  if (f == NULL) {
    fprintf(stderr, "error: can't open '%s' for reading\n", filename);
    abort();
  } 
  char line[1024];
  char seq[10*1024];
  char name[1024];
  int numseqs = 0;
  
  while (fgets(line, 1024, f)) {

    // Strip newline.
    char* s = strchr(line, '\n');
    if (s != NULL) *s = 0;

    if (line[0] == '>'){
      // Header line, i.e. new sequence.
      if (numseqs == 0) {
        // First sequence.
        strcpy(name, line+1);
        seq[0] = 0;
      } else {
        // Save previous sequence.
        add_seq(name, seq);
        seq[0] = 0;
        strcpy(name, line+1);
      }
      numseqs++;      
    } else if (line[0] == '\n') {
      // Empty line, skip
    } else {
      // Assumed sequence.
      if (numseqs == 0){
        fprintf(stderr, "error: invalid FASTA\n");
        fclose(f);
        abort();
      }
      
      // remove any whitespace
      while ((s = strchr(line, ' '))) *s = 0;
      while ((s = strchr(line, '\r'))) *s = 0;
      strcat(seq, line);
    }
  }
  // Save last sequence
  add_seq(name, seq);
  fclose(f);  
  fprintf(stderr, "Read %d sequences.\n", numseqs);
}

static bool is_whitespace(char c)
{
  if (c == ' ' or c == '\n' or c == '\t' or c == '\r') return true;
  return false;
}

static vector<char*> split_string(char* s)
{
  vector<char*> v;
  
  char* t = s;
  bool in_word = false;
  while (*t != 0) {
    if (not is_whitespace(*t)) {
      if (not in_word) v.push_back(t);
      in_word = true;
    } else {
      *t = 0;
      in_word = false;
    }
    t++;
  }
  return v;
}

static inline bool is_gap(char c)
{
  if (c == '-') return true;
  return false;
}

static inline bool is_letter(char c)
{
  char t = toupper(c);
  if (t == 'C' or t == 'G' or t == 'A' or t == 'T' or t == 'N') return true;
  return false;
}

double get_distance(char* seq1, char*seq2, int a, int b)
{
  int n1 = 0, n2 = 0, n = 0;
  for (int i = a; i < b + 1; i++){
    if (is_letter(seq1[i])) n1++;
    if (is_letter(seq2[i])) n2++;
    if (is_letter(seq1[i]) or is_letter(seq2[i]))
      if (toupper(seq1[i]) == toupper(seq2[i])) n++;
  }
  double d = 1.0 - double(n) * 2.0 / double(n1 + n2);
  assert(d >= 0.0);
  return d;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    usage(argv[0]);
    return 1;
  }
  char* filename = argv[1];

  load_fasta(filename);
  int numseqs = seqs.size();

  // Verify all sequences have same length.
  int l = strlen(seqs[0]);
  for (int i = 0; i < numseqs; i++){
    if (int(strlen(seqs[i])) != l) {
      fprintf(stderr,
             "error: sequence '%s' has size %d but sequence '%s' has size %d\n",
              names[i], int(strlen(seqs[i])), names[0], l);
      return 1;
    }
  }

  // Go through sequences and calculate a and b.
  int a = 0, b = l-1;
  for (int i = 0; i < numseqs; i++){
    //fprintf(stderr, "Name: %s\nSeq:  %s\n", names[i], seqs[i]);
    char* seq = seqs[i];
    int t_a = 0, t_b = l - 1;
    // figure out t_a
    while (is_gap(seq[t_a])) t_a++;
    while (is_gap(seq[t_b])) t_b--;
    if (t_a > a) a = t_a;
    if (t_b < b) b = t_b;
  }
  if (a > 10) a = 10;
  if (b < l - 10) b = l - 10;
  fprintf(stderr, "Found a = %d, b = %d\n", a, b);
  
  if (argc == 2) {  
    // Now go through each pair of sequences and print out the distance matrix.
    for (int i = 0; i < numseqs; i++){
      if ((i % 100) == 0 and i != 0){
        fprintf(stderr, "Calculated %d sequence PSI distances\n", i);
      }

      for (int j = i + 1; j < numseqs; j++){
        double d = get_distance(seqs[i], seqs[j], a, b);
        printf("%s %s %.10f\n", names[i], names[j], d); 
      }
    }
  } else if (argc == 3) {
    // If second argument is given (a distance matrix) then print out the
    // distance matrix in normalized form.
    
    // First find average distance in PSI matrix.
    double avg = 0.0; int n = 0;
    for (int i = 0; i < numseqs; i++){
      if ((i % 100) == 0 and i != 0){
        fprintf(stderr, "Calculated %d sequence PSI distances\n", i);
      }
      for (int j = i + 1; j < numseqs; j++){
        avg += get_distance(seqs[i], seqs[j], a, b);
        n++;
      }
    }
    avg = avg / double(n);
    
    fprintf(stderr, "Average PSI distance: %f\n", avg);
    
    // Now find average of the distance matrix.
    FILE* f = fopen(argv[2], "r");
    if (f == NULL) {
      fprintf(stderr, "error: can't open '%s' for reading\n", argv[2]);
    }
    double avgd = 0.0; n = 0;
    char line[1024];
    while (fgets(line, 1024, f)){
      vector<char*> v = split_string(line);
      if (v.size() == 3){
        avgd += atof(v[2]);
        n++;
      } else if (v.size() != 0) {
        fprintf(stderr, "malformed line in '%s'\n", argv[2]);
        abort();
      }
    }
    fclose(f);
    avgd = avgd / double(n);
    fprintf(stderr, "Average matrix distance: %f\n", avgd);
    fprintf(stderr, "Normalization factor: %f\n", avg / avgd);

    fprintf(stderr, "Now normalizing and writing output matrix...\n");
    f = fopen(argv[2], "r");
    if (f == NULL) {
      fprintf(stderr, "error: can't open '%s' for reading\n", argv[2]);
    }
    while (fgets(line, 1024, f)){
      vector<char*> v = split_string(line);
      if (v.size() == 3){
        printf("%s %s %.10f\n", v[0], v[1], atof(v[2]) * avg / avgd);
      }
    }
    fclose(f);
  } else {
    usage(argv[0]);
    return 1;
  }
  return 0;  
}

