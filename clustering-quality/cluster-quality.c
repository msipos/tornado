#include "cfasta.h"

Fasta* fasta;
GArray* clusters;

void load_clustering(const char* filename)
{
  GError* err = NULL;
  GIOChannel* ch = g_io_channel_new_file(filename, "r", &err);
  if (err != NULL){
    g_critical("error while opening clustering file");
    return;
  }
  
  gchar* line;
  while(g_io_channel_read_line(ch, &line, NULL, NULL, NULL) == G_IO_STATUS_NORMAL)
  {
    gchar** cluster_strings = g_strsplit(line, ",", 50000);
    GArray* cluster = g_array_new(FALSE, FALSE, sizeof(Sequence*));
    
    guint n = 0;
    while (cluster_strings[n] != NULL) {
      gchar* seq_id = g_strstrip(cluster_strings[n]);
      Sequence* seq = (Sequence*) g_hash_table_lookup(fasta->hash, 
                                                      seq_id);
      if (seq == NULL){
        g_critical("Sequence '%s' not found!\n", seq_id);
      }
      n++;
      g_array_append_val(cluster, seq);
    }
    //g_message("found %d sequences\n", n);
    
    g_array_append_val(clusters, cluster);
    
    g_free(line);
  }
}

void calculate_CH()
{
  // Using algorithm from Baulik and Bandyopadhyay. See there for
  // notation.

  // First calculate centroid of entire dataset.
  ProbabilitySequence* total_centroid = calc_median(fasta->array);
  //display_prob(total_centroid, 10);
  
  // Calculate TrB - trace of between cluster scatter matrix.
  // k is index over clusters.
  double TrB = 0.0;
  for (guint k = 0; k < clusters->len; k++){
    GArray* cluster = g_array_index(clusters, GArray*, k);
    
    // Calculate centroid of cluster.
    ProbabilitySequence* cluster_centroid = calc_median(cluster);
    //display_prob(cluster_centroid, 10);
    
    double dist = dist_pseq_pseq(cluster_centroid, total_centroid);
    TrB += dist * dist * ((double) cluster->len);
  }
  g_print("TrB = %f\n", TrB);

  // Calculate TrW - trace of within cluster scatter matrix.
  double TrW = 0.0;
  for (guint k = 0; k < clusters->len; k++){
    GArray* cluster = g_array_index(clusters, GArray*, k);

    // Calculate centroid of cluster.
    ProbabilitySequence* cluster_centroid = calc_median(cluster);

    for (guint i = 0; i < cluster->len; i++){
      Sequence* seq = g_array_index(cluster, Sequence*, i);
      double dist = dist_seq_pseq(seq, cluster_centroid);
      TrW += dist*dist;
    }
  }
  g_print("TrW = %f\n", TrW);
  
  // K is number of clusters
  double K = (double) clusters->len;
  
  // n is total num of elements
  double n = (double) fasta->array->len;
  
  double CH = TrB * (n - K) / (K - 1) / TrW;
  g_print("K = %f\n", K);
  g_print("n = %f\n", n);
  g_print("CH = %f\n", CH);
}

int main(int argc, const char** argv)
{
  const char* fasta_filename = argv[1];
  const char* clustering_filename = argv[2];
  
  // Init stuff.
  clusters = g_array_new(FALSE, FALSE, sizeof(GArray*));
  
  GError* err = NULL;
  fasta = load_fasta(fasta_filename, &err);
  if (err != NULL){
    g_critical("error while loading FASTA\n");
    return 1;
  }
  g_print("now loading clustering...\n");
  load_clustering(clustering_filename);
  calculate_CH();
  
  return 0;
}
