#include <LAGraphX.h>
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define MAX_FILENAME 512
#define MAX_LABELS 16

#define RESULTS_DIR "Results/"

#define QUERY_META_FILE "meta.txt"
#define QUERY_SOURCES_FILE "source.txt"
#define QUERY_MAX_STARTING_STATES 16
#define QUERY_MAX_FINAL_STATES 16

#define DO_STRINGIZE(X) #X
#define STRINGIZE(X) DO_STRINGIZE(X)

#define OK(f)                                                             \
  do {                                                                    \
    int res = f;                                                          \
    if (res != 0) {                                                       \
      printf("Error in " #f "() " __FILE__ ":" STRINGIZE(__LINE__) "\n"); \
      printf("LAGraph message %s\n", msg);                                \
      exit(-1);                                                           \
    }                                                                     \
  } while (0)
#define FATAL(...)       \
  do {                   \
    printf(__VA_ARGS__); \
    exit(-1);            \
  } while (0)

#define VERBOSE(...)                  \
  do {                                \
    if (verbose) printf(__VA_ARGS__); \
  } while (0)

void *xmalloc(size_t size) {
  void *res = malloc(size);
  if (res == NULL) exit(-1);
  return res;
}

void *xcalloc(size_t num, size_t size) {
  void *res = calloc(num, size);
  if (res == NULL) exit(-1);
  return res;
}

char msg[LAGRAPH_MSG_LEN];
bool verbose;
GrB_Matrix A;

LAGraph_Graph *GS;

enum query_kind {
  SINGLE_SOURCE,
  SINGLE_DESTINATION,
  ALL_PATHS,
  FIXED_SOURCE_DESTINATION
};

struct query {
  enum query_kind kind;

  int64_t source;
  int64_t dest;

  size_t label_count;
  int64_t labels[MAX_LABELS];
  bool inverse_labels[MAX_LABELS];
  LAGraph_Graph R[MAX_LABELS];

  size_t nqs;
  uint64_t qs[QUERY_MAX_STARTING_STATES];

  size_t nqf;
  uint64_t qf[QUERY_MAX_FINAL_STATES];
};

struct config {
  bool preload;
  bool heatup;
  bool cache_transposed;
  bool profile;

  size_t runs;

  size_t label_count;
  const char *dataset_dir;

  size_t query_count;
  const char *query_dir;
};

struct config DEFAULT_CONFIG = {.preload = true,
                                .heatup = true,
                                .cache_transposed = true,

                                .runs = 5,

                                // Uninitialized values.
                                .label_count = 0,
                                .dataset_dir = NULL,

                                .query_count = 0,
                                .query_dir = NULL};

void init(const struct config *config) {
  // Instead of using conventional LAGraph_Init() use internal
  // LAGr_Init to pass blocking mode.
  // OK(LAGraph_Init(msg));
  OK(LAGr_Init(GrB_BLOCKING, malloc, calloc, realloc, free, msg));

  if (config->profile) GrB_set(GrB_GLOBAL, true, GxB_BURBLE);
}

// Load adjacency matrix for the specified label.
// Returns 0 on success, 1 -- already loaded, -1 -- not loaded.
int load_adjacency_matrix(const struct config *config, size_t label) {
  VERBOSE("Loading adjacency matrix %ld.\n", label);

  assert(config != NULL);

  if (GS[label]) return 1;

  char filename[MAX_FILENAME + 1];

  bool cache_transposed = config->cache_transposed;
  const char *dataset_dir = config->dataset_dir;

  assert(dataset_dir != NULL);

  snprintf(filename, MAX_FILENAME, "%s/%ld.txt", dataset_dir, label);
  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    VERBOSE("Skipped loading adjacency matrix %ld.\n", label);
    return -1;
  }

  OK(LAGraph_MMRead(&A, f, msg));
  GrB_Matrix_wait(A, GrB_MATERIALIZE);
  OK(LAGraph_New(&(GS[label]), &A, LAGraph_ADJACENCY_DIRECTED, msg));

  // A should be moved.
  assert(A == NULL);

  if (cache_transposed) OK(LAGraph_Cached_AT(GS[label], msg));

  fclose(f);

  VERBOSE("Successfully loaded adjacency matrix %ld.\n", label);
  return 0;
}

void load_dataset(struct config *config) {
  VERBOSE("Loading the matrices...\n");

  size_t loaded = 0;

  for (int64_t label = 1; label <= config->label_count; label++) {
    load_adjacency_matrix(config, label);
    if (GS[label] != NULL) loaded++;
  }

  VERBOSE("Successfully loaded %ld adjacency matrices.", loaded);

  size_t total = 0;

  for (int64_t label = 1; label <= config->label_count; label++) {
    if (GS[label]) {
      size_t s = 0;
      GxB_Matrix_memoryUsage(&s, GS[label]->A);
      total += s;

      if (GS[label]->AT != NULL) {
        GxB_Matrix_memoryUsage(&s, GS[label]->AT);
        total += s;
      }
    }
  }

  VERBOSE("Total memory consumption: %ld", total);
  VERBOSE("Loading done!\n");

  VERBOSE("Matrix load has been completed\n");
}

struct query *load_query(const struct config *config, size_t query_number) {
  VERBOSE("Loading query %ld.\n", query_number);

  const char *query_dir = config->query_dir;

  assert(query_dir != NULL);

  GrB_Matrix A = NULL;
  struct query *query = xmalloc(sizeof(struct query));

  char filename[MAX_FILENAME + 1];
  snprintf(filename, MAX_FILENAME, "%s/%ld/" QUERY_META_FILE, query_dir,
           query_number);
  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    VERBOSE("Query %ld is missing. Skipping its loading.\n", query_number);
    return NULL;
  }

  fscanf(f, "%ld", &query->source);
  fscanf(f, "%ld", &query->dest);
  // TODO: Check the source and the destination nodes.

  // Decrease by 1 since MatrixMarket format enumerates rows/cols
  // starting from 0 while GraphBLAS enumerates rows/cols starting
  // from 1.
  query->source--;
  query->dest--;

  // Source or destination equal to -1 (originally 0) means
  // "no fixed source/destination".
  if (query->source == -1 && query->dest == -1)
    query->kind = ALL_PATHS;
  else if (query->dest == -1)
    query->kind = SINGLE_SOURCE;
  else if (query->source == -1)
    query->kind = SINGLE_DESTINATION;
  else
    query->kind = FIXED_SOURCE_DESTINATION;

  fscanf(f, "%ld", &query->nqs);
  assert(query->nqs < QUERY_MAX_STARTING_STATES);

  for (int i = 0; i < query->nqs; i++) {
    fscanf(f, "%ld", &(query->qs[i]));
    query->qs[i]--;
    // TODO: Check starting NFA states.
  }

  fscanf(f, "%ld", &query->nqf);
  assert(query->nqf < QUERY_MAX_FINAL_STATES);
  for (int i = 0; i < query->nqf; i++) {
    fscanf(f, "%ld", &(query->qf[i]));
    query->qf[i]--;
    // TODO: Check final NFA states.
  }

  size_t label_count = 0;
  fscanf(f, "%ld", &label_count);
  for (int i = 0; i < label_count; i++) {
    int64_t label;
    fscanf(f, "%ld", &label);
    // TODO: Check NFA labels.
    query->labels[i] = label;
    query->inverse_labels[i] = label < 0;
  }

  fclose(f);

  // Load NFA adjacency matrices.
  for (int i = 0; i < label_count; i++) {
    snprintf(filename, MAX_FILENAME, "%s/%ld/%ld.txt", query_dir, query_number,
             query->labels[i]);
    FILE *f = fopen(filename, "r");
    if (f == NULL) continue;

    OK(LAGraph_MMRead(&A, f, msg));
    GrB_Matrix_wait(A, GrB_MATERIALIZE);
    OK(LAGraph_New(&(query->R[i]), &A, LAGraph_ADJACENCY_DIRECTED, msg));
    OK(LAGraph_Cached_AT(query->R[i], msg));
    fclose(f);
  }

  VERBOSE("Successfully loaded query %ld.\n", query_number);

  return query;
}

struct query **load_queries(struct config *config) {
  VERBOSE("Loading queries...\n");

  struct query **queries =
      xcalloc(config->query_count + 1, sizeof(struct query));

  size_t loaded = 0;

  // Queries are enumerated starting from 1.
  for (size_t i = 1; i <= config->query_count; i++) {
    queries[i] = load_query(config, i);

    if (queries[i] != NULL) loaded++;
  }

  VERBOSE("Successfully loaded %ld queries\n", loaded);

  return queries;
}

void bench(const struct config *config, struct query **queries) {
  char filename[MAX_FILENAME + 1];

  uint64_t runs = config->runs;
  size_t query_count = config->query_count;
  bool heatup = config->heatup;

  snprintf(filename, MAX_FILENAME, RESULTS_DIR "all.txt");
  FILE *results_f = fopen(filename, "w");
  if (results_f == NULL)
    FATAL("Unable to open all results file by path %s.\n", filename);

  // Zeroth run is a heat-up run.
  for (int run = heatup ? 0 : 1; run <= runs; run++) {
    VERBOSE("Run %d\n", run);
    for (size_t i = 1; i <= query_count; i++) {
      struct query *query = queries[i];

      if (query == NULL) {
        VERBOSE("Query with number %ld isn't present. Skipping.\n", i);
        continue;
      }

      const enum query_kind kind = query->kind;

      LAGraph_Graph G[MAX_LABELS];
      for (size_t i = 0; i <= query->label_count; i++) {
        uint64_t label =
            query->labels[i] > 0 ? query->labels[i] : -query->labels[i];
        int res = load_adjacency_matrix(config, label);
        if (res < 0)
          FATAL("Couldn't find adjacency matrix for label %ld in query %ld\n",
                label, i);

        assert(GS[label] != NULL);
        G[i] = GS[label];
      }

      struct timespec start, finish;

      int64_t S[1] = {query->source};
      int64_t D[1] = {query->dest};

      // Answer vector.
      GrB_Vector reachable = NULL;

      clock_gettime(CLOCK_MONOTONIC, &start);

      switch (kind) {
        case SINGLE_SOURCE:
          OK(LAGraph_2RegularPathQuery(
              &reachable, query->R, query->inverse_labels, query->label_count,
              query->qs, query->nqs, query->qf, query->nqf, G, S, 1, false,
              msg));
          break;
        case SINGLE_DESTINATION:
          OK(LAGraph_2RegularPathQuery(
              &reachable, query->R, query->inverse_labels, query->label_count,
              query->qf, query->nqf, query->qs, query->nqs, G, D, 1, true,
              msg));
          break;
        case ALL_PATHS:
          VERBOSE(
              "Query %ld is ALL PATHS. Such queries aren't supported yet. "
              "Skipping.\n",
              i);
          continue;
      }

      size_t answer = 0;
      GrB_Vector_nvals(&answer, reachable);

      clock_gettime(CLOCK_MONOTONIC, &finish);

      double elapsed = (finish.tv_sec - start.tv_sec) * 1000000.0 +
                       (finish.tv_nsec - start.tv_nsec) / 1000.0;

      printf("%ld,%.0lf,%ld\n", i, elapsed, answer);
      fprintf(results_f, "%ld,%.0lf,%ld\n", i, elapsed, answer);

      OK(GrB_free(&reachable));

      // "If not heat up run..."
      if (run > 0) {
        snprintf(filename, MAX_FILENAME, RESULTS_DIR "%ld.txt", i);
        FILE *f = fopen(filename, "a");
        if (f == NULL) continue;

        fprintf(f, "%.0lf %ld\n", elapsed, answer);
        fclose(f);
      }
    }
  }

  fclose(results_f);
}

void usage() {
  printf(
      "Usage: ./rpq-bench <dataset dir> <label count> <query dir> <query "
      "count>\n");
  printf("\n");
  printf("    -v          Enable verbose logging.\n");
  printf("    -g          Enable performance profiling.\n");
  printf("    -r <runs>   Run count (default: 1).\n");
  printf("    -x          Disable heatup.\n");
  printf("    -p          Disable preloading matricies.\n");
  printf("                Only load required for query evaluation.\n");
  printf("                NB: Also disables total memory consumption.\n");
  printf("    -t          Disable preloading transposed matricies.\n");
  printf("                Enabling the option halves memory consumption\n");
  printf("                but drastically slows 2-RPQ evaluation.\n");
  printf("\n");
  exit(0);
}

void finalize() { LAGraph_Finalize(msg); }

int main(int argc, char **argv) {
  struct config config = DEFAULT_CONFIG;
  opterr = 0;
  char c;

  while ((c = getopt(argc, argv, "xpvgr:t")) != -1) switch (c) {
      case 'p':
        config.preload = false;
        break;
      case 'v':
        verbose = true;
        break;
      case 'g':
        config.profile = true;
        break;
      case 'x':
        config.heatup = false;
        break;
      case 'r':
        config.runs = atoi(optarg);
        break;
      case 't':
        config.cache_transposed = false;
        break;
      case '?':
        if (optopt == 'c')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }

  if (argc - optind < 4) {
    usage();
  }

  config.dataset_dir = argv[optind];
  config.label_count = atoi(argv[optind + 1]);
  config.query_dir = argv[optind + 2];
  config.query_count = atoi(argv[optind + 3]);

  GS = xcalloc(config.label_count + 1, sizeof(LAGraph_Graph));

  VERBOSE("Using dataset dir '%s'\n", config.dataset_dir);
  VERBOSE("Using query dir '%s'\n", config.query_dir);

  init(&config);

  if (config.preload) load_dataset(&config);

  struct query **queries = load_queries(&config);

  bench(&config, queries);

  return 0;
}
