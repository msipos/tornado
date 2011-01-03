#ifndef CFASTA_H
#define CFASTA_H

#include <glib.h>

typedef struct {
  const gchar* seq;
  const gchar* seq_id;
  guint len;
} Sequence;

typedef struct {
  const gchar* seq_id;
  double *data;
  guint len;
} ProbabilitySequence;

typedef struct {
  GHashTable* hash;  // Hash table gchar* seq_id => Sequence*
  GArray* array;
} Fasta;

Fasta* load_fasta(const gchar* filename, GError** err_out);
void display_fasta(Sequence* s);
ProbabilitySequence* calc_median(GArray* arr);
void display_prob(ProbabilitySequence* s, int max);
double dist_seq_pseq(Sequence* s, ProbabilitySequence* ps);
double dist_pseq_pseq(ProbabilitySequence* s1, ProbabilitySequence* s2);
double dist_seq_seq(Sequence* s1, Sequence* s2);

#endif
