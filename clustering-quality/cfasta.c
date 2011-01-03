#include <string.h>
#include "cfasta.h"

void got_sequence(Fasta* fasta, gchar* seq_id, gchar* seq)
{
  //g_message("parsed sequence '%s'", seq_id);
  Sequence* sequence = g_new(Sequence, 1);
  sequence->seq_id = g_intern_string(seq_id);
  sequence->seq = g_intern_static_string(g_ascii_strdown(seq, -1));
  sequence->len = strlen(sequence->seq);
  g_array_append_val(fasta->array, sequence);
  g_hash_table_insert(fasta->hash, (gpointer) sequence->seq_id, sequence);
}

Fasta* load_fasta(const gchar* filename, GError** err_out)
{
  GError* err = NULL;
  GIOChannel* ch = g_io_channel_new_file(filename, "r", &err);
  if (err != NULL){
    *err_out = err;
    return NULL;
  }

  // Initialize Fasta
  Fasta* fasta = g_new(Fasta, 1);
  fasta->array = g_array_new(FALSE, FALSE, sizeof(Sequence*));
  fasta->hash = g_hash_table_new(g_str_hash, g_str_equal);

  gchar* line;
  GString* seq = g_string_new("");
  GString* seq_id = g_string_new("");
  gboolean first = TRUE;
  while(g_io_channel_read_line(ch, &line, NULL, NULL, NULL) == G_IO_STATUS_NORMAL)
  {
    //g_message("line='%s'", line);
    if (line[0] == '>'){
      if (first){
        first = FALSE;
      } else {
        got_sequence(fasta, seq_id->str, seq->str);
      }
      g_string_assign(seq_id, g_strstrip(line + 1));
      g_string_assign(seq, "");
    } else {
      g_string_append(seq, g_strstrip(line));
    }
    g_free(line);
  }
  if (first){
    return NULL;
  } else {
    got_sequence(fasta, seq_id->str, seq->str);
  }

  g_string_free(seq, TRUE);
  g_string_free(seq_id, TRUE);

  g_io_channel_unref(ch);
  return fasta;
}

ProbabilitySequence* calc_median(GArray* arr)
{
  // Verify all sequences have same length.
  guint l = 0;
  for (guint i = 0; i < arr->len; i++){
    Sequence* s = g_array_index(arr, Sequence*, i);
    if (l == 0) l = s->len;
    else {
      if (l != s->len){
        g_critical("Sequences are not padded!");
      }
    }
  }

  ProbabilitySequence* p = g_new(ProbabilitySequence, 1);
  // Construct ProbabilitySequence.
  p->data = g_malloc(l * 6 * sizeof(double));
  p->len = l;

  // Now go through each column.
  for (guint i = 0; i < l; i++){
    guint A=0, C=0, G=0, T=0, g=0, N=0;
    for (guint j = 0; j < arr->len; j++){
      const gchar* seq = g_array_index(arr, Sequence*, j)->seq;
      // Quick'n'dirty
      switch (seq[i]){
        case 'a':
          A++;
          break;
        case 'c':
          C++;
          break;
        case 'g':
          G++;
          break;
        case 't':
          T++;
          break;
        case '-':
          g++;
          break;
        case 'N':
        default:
          N++;
          break;
      }
    }

    double total = (double) arr->len;
    p->data[i*6 + 0] = ((double) A) / total;
    p->data[i*6 + 1] = ((double) C) / total;
    p->data[i*6 + 2] = ((double) G) / total;
    p->data[i*6 + 3] = ((double) T) / total;
    p->data[i*6 + 4] = ((double) g) / total;
    p->data[i*6 + 5] = ((double) N) / total;
  }

  return p;
}

void display_fasta(Sequence* s)
{
  g_print(">%s\n", s->seq_id);
  g_print("%s\n", s->seq);
}

double dist_seq_pseq(Sequence* s, ProbabilitySequence* ps)
{
  g_assert(s->len == ps->len);
  double total = (double) s->len;
  double points = 0.0;
  for (guint i = 0; i < s->len; i++){
    guint j = 0;
    switch(s->seq[i]){
      case 'a':
        j = 0;
        break;
      case 'c':
        j = 1;
        break;
      case 'g':
        j = 2;
        break;
      case 't':
        j = 3;
        break;
      case '-':
        j = 4;
        break;
      case 'n':
      default:
        j = 5;
        break;
    }
    points += ps->data[i*6 + j];
  }
  return 1.0 - points / total;
}

double dist_pseq_pseq(ProbabilitySequence* s1, ProbabilitySequence* s2)
{
  g_assert(s1->len == s2->len);
  double total = (double) s1->len;
  double points = 0.0;
  for (guint i = 0; i < s1->len; i++){
    for (guint j = 0; j < 6; j++){
      double d = (s1->data[i*6 + j] - s2->data[i*6 + j]);
      points += d * d;
    }
  }
  return points / total / 6.0;
}

double dist_seq_seq(Sequence* s1, Sequence* s2)
{
  g_assert(s1->len == s2->len);
  double total = (double) s1->len;
  double points = 0.0;
  for (guint i = 0; i < s1->len; i++){
    if (s1->seq[i] == s2->seq[i]) points += 1.0;
  }
  return 1.0 - points / total;
}

void display_prob(ProbabilitySequence* s, int max)
{
  for (guint j = 0; j < 6; j++){
    gchar c;
    switch(j){
      case 0:
        c = 'A';
        break;
      case 1:
        c = 'C';
        break;
      case 2:
        c = 'G';
        break;
      case 3:
        c = 'T';
        break;
      case 4:
        c = '-';
        break;
      case 5:
        c = 'N';
        break;
    }
    g_print("%c:", c);
    for (guint i = 0; i < s->len; i++){
      g_print(" %.2f", s->data[i * 6 + j]);
      if (i > max) break;
    }
    g_print("\n");
  }
}
