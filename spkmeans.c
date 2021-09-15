#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include "spkmeans.h"

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    printf("Invalid Input!");
  }
  k = atoi(argv[1]);              /* get k from cmd*/
  interface(k, argv[2], argv[3]); /* goal and file from cmd*/
  return 0;
}

void interface(int k, char *goal, char *filename)
{
  int *d_n = get_d_n(filename); /* obtain N and the dimension of
   the points in the file*/
  n = d_n[1];
  if (k >= n)
  {
    exit_on_input();
  }
  d = d_n[0];
  data = build_arr(d_n[0], d_n[1], filename); /* build datapoints
   from the input file*/
  W = build_w(data, n);
  D = build_D(W, n);
  free_point(data, n);
  if (strcmp(goal, "wam") == 0) /* print W*/
  {
    print_and_free_mat(W, n);
    free_mat(D, n);
  }
  if (strcmp(goal, "ddg") == 0) /* Print the diagonal*/
  {
    print_and_free_mat(D, n);
    free_mat(W, n);
  }
  if (strcmp(goal, "lnorm") == 0) /* Print Laplacian*/
  {
    L = build_L(W, D, n);
    print_and_free_mat(L, n);
    free_mat(D, n);
    free_mat(W, n);
  }

  if (strcmp(goal, "jacobi") == 0) /* Calculate jacobi
   and print eigenvectors and eigengaps*/
  {
    arr_d = build__double_arr(n, filename);
    /*tup=L,V*/
    tup = jacobi(arr_d, n);
    eigengap = get_eigengap(tup[0], n);
    for (i = 0; i < n - 1; i++)
    {
      printf("%.4f,", round(eigengap[i].eigen));
    }
    printf("%.4f\n", round(eigengap[n - 1].eigen));
    print_jacobi(tup[1], n);
    free(eigengap);
    free_mat(D, n);
    free_mat(W, n);
    free_mat(tup[0], n);
    free(tup);
  }

  if (strcmp(goal, "spk") == 0) /* Perform the full spectral clustring
   and Kmeans and output the final centroids*/
  {
    run_spk(W, D, n, k);
  }
}

/*Exit and terminate after memory error*/
void exit_run()
{
  printf("An Error Has Occured");
  exit(-1);
}
/*Exit and terminate when the input is invalid*/
void exit_on_input()
{
  printf("Invalid Input!");
  exit(-1);
}

/*Print matrix from algorithm and free it*/
void print_and_free_mat(double **mat, int n)
{
  for (i = 0; i < n; i++)
  {
    for (q = 0; q < n - 1; q++)
    {
      printf("%.4f", round(mat[i][q]));
      printf(",");
    }
    if (i==n-1){
      printf("%.4f", round(mat[i][n - 1]));
    }
    else {
    printf("%.4f\n", round(mat[i][n - 1])); }
    free(mat[i]);
  }
  free(mat);
}
/*Run the full spk algorithm*/
void run_spk(double **W, double **D, int n, int k)
{
  L = build_L(W, D, n); /*Laplacian*/
  tup = jacobi(L, n);
  L = tup[0];
  V = tup[1];
  free_mat(W, n);
  free_mat(D, n);
  eigengap = get_eigengap(L, n);/*Get gigengaps from L diagonal*/
  qsort(eigengap, n, sizeof(Eigengap), compare_eigen); /*sort the eigengaps
   and keep the indices for later use*/
  free_mat(L, n);
  if (k == 0)
  {
    k = get_k(eigengap, n); /* Run the Eigengep Heuristic*/
  }

  U = sort_U(V, eigengap, n, k); /*the matrix U*/
  print_jacobi(V,n);
  printf("----\n");
  T = getT(U, n, k);             /* Get the matrix T as the datapoints
   - K means will be performed on T*/
  free_mat(U, n);
  free_mat(V, n);
  free(eigengap);
  free(tup);
  datapoints = build_DataPoints(T, n, k);
  clusters = create_clusters(datapoints, k, k); /* first k datapoints
  are the initial centroids*/
  k_means(datapoints, clusters, k, k, n, 300); /*Run the Kmeans algorithm on 
  the datapoints and print the result*/
  for (i = 0; i < n; i++)
  {
    free(datapoints[i].coords);
  }
  free(datapoints);
  free(clusters);
  free_mat(T, n);
}

/* Print  the jacobi output and free matrix*/
void print_jacobi(double **arr, int n)
{

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n - 1; j++)
    {
      printf("%.4f", round(arr[j][i]));
      printf(",");
    }
    if (i==n-1){
      printf("%.4f", round(arr[n - 1][i]));
    }
    else{
      printf("%.4f\n", round(arr[n - 1][i]));
    }
  }
  free_mat(arr, n);
}
/*Free double** NxN matrix */
void free_mat(double **arr, int n)
{

  for (i = 0; i < n; i++)
  {
    free(arr[i]);
  }
  free(arr);
}
/*Free point* NxN matrix */
void free_point(point **arr, int n)
{

  for (i = 0; i < n; i++)
  {
    free(arr[i]->coords);
    free(arr[i]);
  }
  free(arr);
}

/*Build the Wrighted adjacency Matrix*/
double **build_w(point **points, int n)
{
  d = points[0]->dim;
  W = (double **)malloc(sizeof(double *) * n);
  if (W == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    W[i] = (double *)malloc(sizeof(double) * n);
    if (W[i] == NULL)
    {
      exit_run();
    }
  }
  for (i = 0; i < n; i++)
  {
    point_i = points[i];
    for (j = 0; j < n; j++)
    {

      weight = 0.0;
      if (i == j)
      {
        W[i][i] = 0.0;
      }
      else
      {
        point_j = points[j];
        for (b = 0; b < d; b++)
        {
          weight += pow(point_i->coords[b] - point_j->coords[b], 2);
        }
        weight = pow(weight, 0.5);
        weight = weight * -0.5;
        weight = exp(weight);
        W[i][j] = weight;
      }
    }
  }
  return W;
}

/*Build the Diagonal Degree Matrix*/
double **build_D(double **W, int n)
{
  D = (double **)malloc(sizeof(double *) * n);
  if (D == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    D[i] = (double *)calloc(n, sizeof(double));
    if (D[i] == NULL)
    {
      exit_run();
    }
  }
  for (i = 0; i < n; i++)
  {
    sum = 0.0;
    for (j = 0; j < n; j++)
    {
      sum += W[i][j];
    }
    D[i][i] = sum;
  }

  return D;
}
/*Build the Diagonal Degree Matrix ^-0.5*/
double **get_D_0_5(double **D, int n)
{

  D_05 = (double **)malloc(sizeof(double *) * n);
  if (D_05 == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    D_05[i] = (double *)calloc(n, sizeof(double));
    if (D_05[i] == NULL)
    {
      exit_run();
    }
  }
  for (i = 0; i < n; i++)
  {
    D_05[i][i] = 1 / pow(D[i][i], 0.5);
  }

  return D_05;
}
/*Get the laplacian from W and D*/
double **build_L(double **W, double **D, int n)
{
  double **L;
  /* Lnorm = I-D^-0.5*W^D-0.5: */
  unit = get_unit(n); /*Unit matrix*/
  d_05 = get_D_0_5(D, n);
  temp1 = mult_diag(d_05, W, n, 1);
  temp2 = mult_diag(temp1, d_05, n, 2);
  L = min_mat(unit, temp2, n);
  free_mat(unit, n);
  free_mat(temp1, n);
  free_mat(temp2, n);
  free_mat(d_05, n);
  return L;
}

/* multiply two nXn matrixes, one or two of them are diagonal*/
double **mult_diag(double **mat1, double **mat2, int n, int diag)
{ 
/*diag indicates which one is the diagonal, 1 is mat1,
 2 is 2mat2, and 3 is both*/
  result = (double **)malloc(sizeof(double *) * n);
  if (result == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    result[i] = (double *)calloc(n, sizeof(double));
    if (result[i] == NULL)
    {
      exit_run();
    }
  }

  if (diag == 3)
  { /* both are diagonal*/

    for (i = 0; i < n; i++)
    {
      result[i][i] = mat1[i][i] * mat2[i][i];
    }
  }

  else if (diag == 1)
  { /* the first mat is diagonal*/
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        result[i][j] = mat1[i][i] * mat2[i][j];
      }
    }
  }
  else
  { /*the second mat is diagonal*/
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        result[i][j] = mat1[i][j] * mat2[j][j];
      }
    }
  }

  return result;
}

/* return mat1-mat2*/
double **min_mat(double **mat1, double **mat2, int n)
{
  result = (double **)malloc(sizeof(double *) * n);
  if (result == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    result[i] = (double *)calloc(n, sizeof(double));
    if (result[i] == NULL)
    {
      exit_run();
    }
  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      result[i][j] = mat1[i][j] - mat2[i][j];
    }
  }

  return result;
}
/* Get the identity matrix*/
double **get_unit(int n)
{
  result = (double **)malloc(sizeof(double *) * n);
  if (result == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    result[i] = (double *)calloc(n, sizeof(double));
    assert(result[i] != NULL);
    if (result[i] == NULL)
    {
      exit_run();
    }
  }

  for (i = 0; i < n; i++)
  {
    result[i][i] = 1.0;
  }

  return result;
}
/* return dimension(d) and size (N) of datapoints*/
int *get_d_n(char *filename)
{
  count = 0;
  d = 0;
  count_d = 1;
  fp = fopen(filename, "r");
  if (fp == NULL)
  {
    exit_on_input();
  }
  while (fscanf(fp, "%lf%c", &value, &c) == 2)
  {
    if (count_d)
    {
      d++;
    }
    if (c == '\n' || c == '\r')
    {
      count++;
      if (count_d)
      {
        count_d = 0;
      }
    }
  }
  if (c != '\n' && c != '\r')
  {
    count++;
  }
  fclose(fp);
  n = count;
  r[0] = d;
  r[1] = n;
  return r;
}
/* build data_points as point*struct*/
point **build_arr(int d, int n, char *filename)
{
  arr = (point **)malloc(sizeof(point *) * n);
  if (arr == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    arr[i] = (point *)malloc(sizeof(point));
    if (arr[i] == NULL)
    {
      exit_run();
    }
  }
  fp = fopen(filename, "r");
  if (fp == NULL)
  {
    exit_on_input();
  }
  for (i = 0; i < n; i++)
  {
    arr[i]->dim = d;
    arr[i]->coords = (double *)malloc(d * sizeof(double));
    if (arr[i]->coords == NULL)
    {
      exit_run();
    }
  }

  j = 0;
  count = 0;
  while (fscanf(fp, "%lf%c", &value, &c) == 2)
  {

    arr[count]->coords[j] = value;

    j++;
    if (c == '\n' || c == '\r')
    {
      count++;
      j = 0;
    }
  }
  if (c != '\n' && c != '\r')
  {
    arr[count]->coords[j] = value;
  }

  fclose(fp);

  return arr;
}
/* build data_points as double ** arr only for symetric input for Jacobi */
double **build__double_arr(int n, char *filename)
{
  double **arr_d = (double **)malloc(sizeof(double *) * n);
  if (arr_d == NULL)
  {
    exit_run();
  }

  for (i = 0; i < n; i++)
  {
    arr_d[i] = (double *)malloc(sizeof(double) * n);
    if (arr_d[i] == NULL)
    {
      exit_run();
    }
  }
  fp = fopen(filename, "r");
  if (fp == NULL)
  {
    exit_on_input();
  }
  j = 0;
  count = 0;

  while (fscanf(fp, "%lf%c", &value, &c) == 2)
  {

    arr_d[count][j] = value;

    j++;
    if (c == '\n' || c == '\r')
    {
      count++;
      j = 0;
    }
  }
  if (c != '\n' && c != '\r')
  {
    arr_d[count][j] = value;
  }
  fclose(fp);

  return arr_d;
}

/* Pivot inside the jacobi part according to the algorithm*/
double *pivot(double **arr, int n)
{
  max = 0;
  best_i = 0;
  best_j = 0;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (i != j && fabs(arr[i][j]) > max)
      {
        best_i = i;
        best_j = j;
        max = fabs(arr[i][j]);
      }
    }
  }

  if (best_j == 0 && best_i == 0)
  {
    result3[0] = -1;
  }
  else
  {
    result3[0] = max;
  }
  result3[1] = best_i;
  result3[2] = best_j;
  return result3;
}
/*Get C and S for jacobi*/
double *get_c_s(double **arr, int i, int j)
{
  theta = (arr[j][j] - arr[i][i]) / (2 * arr[i][j]);
  t = sign(theta) / (fabs(theta) + sqrt(theta * theta + 1));
  c_2 = 1 / (sqrt(t * t + 1)); /*C*/
  s = t * c_2;
  result2[0] = c_2;
  result2[1] = s;
  return result2;
}

/*Get the sign of a double*/
int sign(double x)
{
  if (x >= 0)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}

/*In-place multiplication of rotation matrix P for jacobi*/
void mult_p(double **mat1, int i, int j, int n, double c, double s)
{
  for (x = 0; x < n; x++)
  {
    double temp_i = mat1[x][i];
    double temp_j = mat1[x][j];
    mat1[x][i] = c * temp_i - s * temp_j;
    mat1[x][j] = s * temp_i + c * temp_j;
  }
}

/*Get P for jacobi*/
double **get_P(double c, double s, int i, int j, int n)
{
  p = get_unit(n);
  p[i][i] = c;
  p[j][j] = c;
  p[j][i] = -s;
  p[i][j] = s;
  return p;
}

/* step from A to A' - jacobi */
double **get_a_tag(double **arr, int n, double c, double s, int i, int j)
{
  result = (double **)malloc(sizeof(double *) * n);
  if (result == NULL)
  {
    exit_run();
  }
  for (d = 0; d < n; d++)
  {
    result[d] = (double *)malloc(n * sizeof(double));
    if (result[d] == NULL)
    {
      exit_run();
    }
  }

  for (q = 0; q < n; q++)
  {
    if (q != i && q != j)
    {
      result[q][i] = c * arr[q][i] - s * arr[q][j];
      result[q][j] = (c * arr[q][j] + s * arr[q][i]);
      result[i][q] = result[q][i];
      result[j][q] = result[q][j];
    }
  }

  result[i][i] =
      (c * c * arr[i][i]) + (s * s * arr[j][j]) - (2 * s * c * arr[i][j]);

  result[j][j] =
      pow(s, 2) * arr[i][i] + pow(c, 2) * arr[j][j] + 2 * s * c * arr[i][j];
  result[i][j] = 0;
  result[j][i] = 0;

  for (d = 0; d < n; d++)
  {
    for (f = 0; f < n; f++)
    {
      if (d != i && d != j && f != i && f != j)
      {
        result[d][f] = arr[d][f];
      }
    }
  }

  return result;
}
/* forbenius norm for jacobi*/
double forbenius_norm(double **arr, int n)
{
  target = 0.0;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      target += (arr[i][j] * arr[i][j]);
    }
  }
  return target;
}
/*Get the square of the diag - Jacobi*/
double square_diag(double **arr, int n)
{
  sq = 0.0;
  for (i = 0; i < n; i++)
  {
    sq += (arr[i][i] * arr[i][i]);
  }
  return sq;
}

/* Calculate "off" for jacobi */
double off(double **arr, int n)
{
  return forbenius_norm(arr, n) - square_diag(arr, n);
}
/*  Check convergence of jacobi algorithm*/
int check_convergence(double x, double y)
{
  if (x - y <= 0.000000000000001)
  {
    return 1;
  }
  return 0;
}

/* The full jacobi algorithm */
double ***jacobi(double **L, int n)
{
  double **V;
  first = 1; /*Indicator for the first iteration*/
  count = 0; /*Iterations counter*/
  converged = 0; /*Indicator for the convergence*/
  V = L; /*Only to avoid python warnings*/
  while (!converged && count < 100)
  {
    piv = pivot(L, n);
    best_i = piv[1];
    best_j = piv[2];
    if (piv[0] == -1)
    {
      converged = 1;
      L_2 = L;
      break;
    }
    c_s = get_c_s(L, best_i, best_j); /*Get C and S*/

    if (first) /* this is the first iteration - V is P*/
    {
      P = get_P(c_s[0], c_s[1], best_i, best_j, n);
      V = P;
      first = 0;
    }
    else
    {
      mult_p(V, best_i, best_j, n, c_s[0], c_s[1]); /*Update V by multiplying P with the current V*/
    }
    count++;
    of_a = off(L, n);
    temp3 = L;
    L = get_a_tag(L, n, c_s[0], c_s[1], best_i, best_j); /*Step from A to A'*/
    free_mat(temp3, n);
    of_a_tag = off(L, n);
    if (check_convergence(of_a, of_a_tag))
    {
      L_2 = L;
      converged = 1;
    }
    if (count == 100)
    {
      L_2 = L;
    }
  }
  result_full = (double ***)malloc(sizeof(double **) * 2);
  if (result_full == NULL)
  {
    exit_run();
  }
  result_full[0] = L_2;
  result_full[1] = V;
  return result_full; /*Teturn Laplacian and the eigenvectors*/
}

/*Get the eigangaps from the diagonal*/
Eigengap *get_eigengap(double **arr, int n)
{
  eigengap = (Eigengap *)malloc(sizeof(Eigengap) * n);
  if (eigengap == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    eigengap[i].eigen = arr[i][i];
    eigengap[i].index = i; /* save the index in order 
  to get it later after sort*/
  }
  return eigengap;
}

/* Compare two eigengaps structs with indexes */
int compare_eigen(const void *a, const void *b)
{
  if ((*(Eigengap *)a).eigen > (*(Eigengap *)b).eigen)
    return 1;
  else if ((*(Eigengap *)a).eigen < (*(Eigengap *)b).eigen)
    return -1;
  else
    return 0;
}

/*The k method Heristic according to the algorithm*/
int get_k(Eigengap *arr, int n)
{
  arg_max = 0;
  i_max = 0;
  for (i = 0; i < n / 2; i++)
  {
    curr = fabs(arr[i].eigen - arr[i + 1].eigen);
    if (curr > arg_max) /*Update max*/
    {
      arg_max = curr;
      i_max = i;
    }
  }
  return i_max + 1; /* counting from 1*/
}

/*Sort the matrix V according to the eigengaps order using their indices*/
double **sort_U(double **V, Eigengap *eigens, int n, int k)
{
  result = (double **)malloc(sizeof(double *) * n);
  if (result == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    result[i] = (double *)malloc(sizeof(double) * k);
    if (result[i] == NULL)
    {
      exit_run();
    }
  }
  for (i = 0; i < k; i++)
  {
    int ind = eigens[i].index;
    for (j = 0; j < n; j++)
    {
      result[j][i] = V[j][ind];
    }
  }
  return result;
}

/* Get T by normalising U*/
double **getT(double **U, int n, int k)
{
  T = (double **)malloc(sizeof(double *) * n);
  if (T == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    T[i] = (double *)malloc(k * sizeof(double));
    if (T[i] == NULL)
    {
      exit_run();
    }
  }
  for (i = 0; i < n; i++)
  {
    sum = 0.0;
    for (b = 0; b < k; b++)
    {
      sum += pow(U[i][b], 2);
    }
    if (sum >= 0)
    {
      sum = pow(sum, 0.5);
    }
    else /*Only for an edge case in which the zero vector is created*/
    {
      sum = sum;
    }
    for (j = 0; j < k; j++)
    {
      if (sum == 0)/*Only for an edge case in which the zero vector is
       created*/
      {
        T[i][j] = U[i][j];
      }
      else
      {
        T[i][j] = U[i][j] / sum;
      }
    }
  }
  return T;
}

/* Get the data points as point struct, for k means*/
point *build_DataPoints(double **T, int n, int k)
{

  datapoints = (point *)malloc(sizeof(point) * n);
  if (datapoints == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    datapoints[i].dim = k;

    datapoints[i].coords = (double *)malloc(sizeof(double) * k);

    if (datapoints[i].coords == NULL)
    {
      exit_run();
    }
  }
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < k; j++)
    {
      datapoints[i].coords[j] = T[i][j];
    }
  }

  return datapoints;
}

/* Create initial first k centroids */
cluster *create_clusters(point *data_points, int k, int d)
{
  clusters = (cluster *)malloc(sizeof(cluster) * k);
  if (clusters == NULL)
  {
    exit_run();
  }
  /* initialize the clusters:*/
  for (i = 0; i < k; i++)
  {
    sum_arr = (double *)malloc(sizeof(double) * d);
    if (sum_arr == NULL)
    {
      exit_run();
    }

    for (b = 0; b < d; b++)
    {
      sum_arr[b] = 0.0;
    }

    clusters[i].cent = data_points[i];
    clusters[i].sum = sum_arr;
    clusters[i].size = 0.0;
  };

  return clusters;
}

/*K means algorithm fromm HW1 and HW2*/
void k_means(point *data_points, cluster *clusters, int k, int d, int n,
             int max_iter)
{
  int iter_count = 0;
  for (i = 0; i < max_iter; i++)
  { /*max_iter*/
    for (j = 0; j < n; j++)
    { /* number of points*/
      point point_x;
      point_x = data_points[j];
      arg_min = 0.0;
      for (b = 0; b < d; b++)
      {
        arg_min = arg_min + (point_x.coords[b] - clusters[0].cent.coords[b]) *
        (point_x.coords[b] - clusters[0].cent.coords[b]);
      }
      assign_to = 0;
      for (y = 0; y < k; y++)
      {
        point centroid = clusters[y].cent;
        double *coords_cent;
        coords_cent = centroid.coords;
        euclid_dist = 0.0;
        for (b = 0; b < d; b++)
        {
          euclid_dist = euclid_dist + (point_x.coords[b] - coords_cent[b]) *
           (point_x.coords[b] - coords_cent[b]);
        }

        if (euclid_dist < arg_min)
        {

          arg_min = euclid_dist;
          assign_to = y;
        }
      }
      /* assign here*/
      coords_x = point_x.coords;
      for (q = 0; q < d; q++)
      {
        clusters[assign_to].sum[q] += coords_x[q];
      }
      clusters[assign_to].size += 1.0;
    }
    changed = 0;
    for (x = 0; x < k; x++)
    {

      sums_and_divide = (double *)malloc(d * sizeof(double));
      if (sums_and_divide == NULL)
      {
        exit_run();
      }
      for (b = 0; b < d; b++)
      {
        if (clusters[x].size > 0)
        {
          sums_and_divide[b] = clusters[x].sum[b] / clusters[x].size;
        }
        else
        {
          sums_and_divide[b] = clusters[x].cent.coords[b];  /*Only for an edge
           case in which a cluster remained empty*/
        }
      }

      curr_cent_node = clusters[x].cent.coords;
      for (z = 0; z < d; z++)
      {

        if (curr_cent_node[z] != sums_and_divide[z])
        {
          changed = 1;
          break;
        }
      }
      if (i != 0)
      {
        free(clusters[x].cent.coords);
      }
      clusters[x].cent.coords = sums_and_divide;
    }
    iter_count++;

    if (changed == 0)

    {
      break;
    }
    for (q = 0; q < k; q++)
    {

      sum_arr = (double *)malloc(sizeof(double) * d);
      if (sum_arr == NULL)
      {
        exit_run();
      }

      for (b = 0; b < d; b++)
      {
        sum_arr[b] = 0.0;
      }
      free(clusters[q].sum);
      clusters[q].sum = sum_arr;
      clusters[q].size = 0.0;
    }
  }
  for (i = 0; i < k; i++)
  {
    free(clusters[i].sum);
    cent = clusters[i].cent.coords;
    for (q = 0; q < d - 1; q++)
    {
      printf("%.4f", round(cent[q]));
      printf(",");
    }
    if (i==k-1){
          printf("%.4f", round(cent[d - 1]));
    }
    else{
      printf("%.4f\n", round(cent[d - 1]));
    }
    free(clusters[i].cent.coords);
  }
}