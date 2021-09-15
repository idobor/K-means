double **result;
double result2[2];
double result3[3];
int r[2];

#ifndef SPKMEANS_H_
#define SPKMEANS_H_
#define round(x) ((((x) > (-0.00005)) && ((x) < (0))) ? (0.0) : (x))

typedef struct
{
  double *coords;
  int dim;
} point;

typedef struct
{

  double eigen;
  int index;
} Eigengap;

typedef struct
{
  point cent;
  double *sum;
  float size;
} cluster;

double **build_w(point **points, int n);

double **build_D(double **W, int n),
    **get_D_0_5(double **D, int n), **mult_diag(double **mat1, double **mat2, int n, int diag), **build_L(double **W, double **D, int n),
    **get_unit(int n), *pivot(double **arr, int n),
    *get_c_s(double **arr, int i, int j),
    **get_P(double c, double s, int i, int j, int n),
    **get_a_tag(double **arr, int n, double c, double s, int i, int j),
    square_diag(double **arr, int n),
    ***jacobi(double **L, int n), **getT(double **U, int n, int k),
    **sort_U(double **V, Eigengap *eigens, int n, int k),
    **min_mat(double **mat1, double **mat2, int n),
    forbenius_norm(double **arr, int n), off(double **arr, int n),
    **build__double_arr(int n, char *filename);

double arg_max, value, **W, **U, **D, **L, **D_05, ***result_full, weight,
    **unit, **d_05, **arr_d, **temp1, **temp2, temp, **p, **P, *piv, **T, *sum_arr,
    arg_min, euclid_dist, *coords_x, *sums_and_divide, *curr_cent_node, max, target, sq,
    *cent, *c_s, sum, theta, t, s, c_2, **V, ***tup, **L_2, of_a, of_a_tag, **temp3, curr;

int check_convergence(double x, double y),
    compare_eigen(const void *a, const void *b), get_k(Eigengap *arr, int n),
    sign(double x), *get_d_n(char *filename);

int count, d, n, q, b, count_d, i_max, assign_to, changed, best_i, best_j, x, y,
    first, count, converged, i, j, d, z, k, f, first_kmeans_iter;

Eigengap *get_eigengap(double **arr, int n), *eigengap;
point *build_DataPoints(double **T, int n, int k), *datapoints,
    **build_arr(int d, int n, char *filename), **arr, **data, *point_i, *point_j;

cluster *create_clusters(point *data_points, int k, int d), *clusters;

void k_means(point *data_points, cluster *clusters, int k, int d, int n,
             int max_iter);
void print_jacobi(double **arr, int n);
void free_mat(double **arr, int n);
void free_point(point **arr, int n);
void mult_p(double **mat1, int i, int j, int n, double c, double s);
void print_and_free_mat(double **mat, int n);
void run_spk(double **W, double **D, int n, int k);
void exit_run();
void exit_on_input();

void interface(int k, char *goal, char *filename);

char c;

FILE *fp;

#endif