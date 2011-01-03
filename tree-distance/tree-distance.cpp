#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <map>
using std::vector;

struct TreeNode {
  char* name;
  double depth;
  vector<TreeNode*> children;
  bool marked;
  TreeNode* parent;
};

//////////////////////////////////////////////////////////////////////
// LOADING A NEWICK FILE
//////////////////////////////////////////////////////////////////////

double read_double(char** s)
{
  char buf[1024];
  char* b = buf;
  while (isdigit(**s) or **s == '.' or **s == 'E' or **s == 'e'
         or **s == '-'){
    *b = **s;
    b++; (*s)++;
  }
  *b = 0;

  return atof(buf);
}

TreeNode* load_node(char** s)
{
  TreeNode* node = new TreeNode;
  node->name = NULL;
  node->depth = 0.0;
  node->marked = false;
  node->parent = NULL;

  if (**s == '(') {
    // Not an leaf node.
    (*s)++;
    for (;;){
      TreeNode* child = load_node(s);
      node->children.push_back(child);
      child->parent = node;
      if (**s == ','){
        (*s)++;
        continue;
      } else if (**s == ')') {
        (*s)++;
        char buf[1024];
        char* b = buf;
        while (**s != ':' and **s != 0){
          *b = **s;
          b++;
          (*s)++;
        }
        *b = 0;
        node->name = strdup(buf);
        break;
      } else {
        assert(true);
      }
    }
  } else {
    char buf[1024];
    char* b = buf;
    while (**s != ':' and **s != ',' and **s != ')'){
      *b = **s;
      b++; (*s)++;
    }
    *b = 0;
    node->name = strdup(buf);
    // printf("read name: %s\n", buf);
  }    
  if (**s == ':'){
    (*s)++;
    node->depth = read_double(s);
  }
  return node;
}

TreeNode* load_newick(FILE* f)
{
  fseek(f, 0, SEEK_END);
  size_t len = (size_t) ftell(f);
  fseek(f, 0, SEEK_SET);
  char* buffer = (char*) malloc(len+1);
  size_t rd = fread(buffer, 1, len, f);
  if (rd != len) {
    fprintf(stderr, "error: read error while loading newick file.\n");
    return NULL;
  }
  buffer[len] = 0;
  
  char* s = buffer;
  return load_node(&s);  
}

void display_newick(TreeNode* node, int offset)
{
  for (int i = 0; i < offset; i++){
    printf(" ");
  }
  if (node->name != NULL) printf("+-%s, %f\n", node->name, node->depth);
  else printf("+-NULL\n");
  if (node->children.size() > 0){
    for (int i = 0; i < offset+1; i++){
      printf(" ");
    }
    printf("|\n");
    for (unsigned int i = 0; i < node->children.size(); i++){
      display_newick(node->children[i], offset+1);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// TREE OPERATIONS
//////////////////////////////////////////////////////////////////////

int count_leaves(TreeNode* node)
{
  if (node->children.size() == 0)
    return 1;

  int l = 0;
  for (unsigned int i = 0; i < node->children.size(); i++){
    l+= count_leaves(node->children[i]);
  }
  return l;
}

double max_depth(TreeNode* node)
{
  double d = 0.0;
  for (unsigned int i = 0; i < node->children.size(); i++){
    double dc = max_depth(node->children[i]);
    if (dc > d) d = dc;
  }
  return d + node->depth;
}

// Calculate all distances coming from this node.
void calc_distances_r(TreeNode* center, TreeNode* node, double distance, TreeNode* came_from)
{
  // Node is not marked and it's a leaf node
  if ((not node->marked) and (node->children.size() == 0)) 
    printf("%s %s %.10f\n", center->name, node->name, distance);
  
  // Recurse into children
  for (int i = 0; i < int(node->children.size()); i++){
    TreeNode* child = node->children[i];
    if (came_from != child) {
      calc_distances_r(center, child, distance + child->depth, node);
    }
  }
  
  // Recurse into parents
  TreeNode* parent = node->parent;
  if ((came_from != parent) and (parent != NULL)) {
    calc_distances_r(center, parent, distance + node->depth, node);
  }
}

void do_distance_r(TreeNode* n)
{
  n->marked = true;
  if (n->children.size() == 0){
    calc_distances_r(n, n, 0.0, NULL);
  }
  
  for (int i = 0; i < int(n->children.size()); i++){
    do_distance_r(n->children[i]);
  }
}


//////////////////////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////////////////////

int main(int argc, char**argv)
{
  if (argc < 2){
    fprintf(stderr, "error: need newick file.\n");
  }

  // Load the newick file
  char* filename = argv[1];
  FILE* f = fopen(filename, "r");
  TreeNode* root = load_newick(f);
  fclose(f);
  
  //display_newick(root, 0);
  
  fprintf(stderr, "number of leaves: %d\n", count_leaves(root));
  
  do_distance_r(root);
}
