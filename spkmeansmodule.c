#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <Python.h>
#include "spkmeans.h"
PyObject *k_py;
PyObject *listObj;
int i, j, d, b;
double **W;
double **D;
double **L;
double **D_05;
double ***tup;
double **T;
cluster *clusters;
point *datapoints;
point **data;
double **V;
double *sum_arr;
static PyObject *T_Python(int n, double **T, int k);
cluster *create_py_clusters(point *data_points, int k, int d, double *keys);
static PyObject *standard_interface(PyObject *self, PyObject *args);
static PyObject *spk_a(PyObject *self, PyObject *args);
static PyObject *spk_b(PyObject *self, PyObject *args);
double **T_from_python(PyObject *dataObj, int n, int k);
static PyObject *spk_python_1(int k, char *filename);
void spk_python_2(double **T, double *keys, int n, int k);
double *get_keys_from_py(PyObject *listObj, int k);

/*interface for C-API python, for non spk goals*/
static PyObject *standard_interface(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &listObj))
  {
    printf("An Error Has Occured");
    return NULL;
  }
  k_py = PyList_GetItem(listObj, 0); /*k*/
  long k_long = PyLong_AsLong(k_py); /*Convert python k to C k*/
  int k = (int)k_long;
  char *goal;
  if (!PyArg_Parse(PyList_GetItem(listObj, 1), "s", &goal))
  {
    printf("An Error Has Occured");
    return NULL;
  }
  char *filename;
  if (!PyArg_Parse(PyList_GetItem(listObj, 2), "s", &filename))
  {
    printf("An Error Has Occured");
    return NULL;
  }
  interface(k, goal, filename); /*regular C interface*/
  return Py_BuildValue("");     /*No return value*/
}

/*First C-API python inteface only for spk - argumnets from cmd*/
static PyObject *spk_a(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &listObj))
  {
    printf("An Error Has Occured");
    return NULL;
  }
  PyObject *k_py = PyList_GetItem(listObj, 0);
  long k_long = PyLong_AsLong(k_py); /*Convert python k to C k*/
  int k = (int)k_long;
  char *goal;
  if (!PyArg_Parse(PyList_GetItem(listObj, 1), "s", &goal))
  {
    printf("An Error Has Occured");
    return NULL;
  }
  char *filename;
  if (!PyArg_Parse(PyList_GetItem(listObj, 2), "s", &filename))
  {
    printf("An Error Has Occured");
    return NULL;
  }
  return spk_python_1(k, filename);
}
/*Second C-API python inteface only for spk -
 argumnets are T and the keys selected by kmeans++ from python*/
static PyObject *spk_b(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &listObj))
  {
    printf("An Error Has Occured");
    return NULL;
  }
  /* keys are the first element, T is the second*/
  int k = PyList_Size(PyList_GetItem(listObj, 0));
  PyObject *dataObj = PyList_GetItem(listObj, 1); /*The datapoints from python*/
  int n = PyList_Size(dataObj);
  double *keys = get_keys_from_py(PyList_GetItem(listObj, 0), k); /*Convert
   the keys from Pyobject to C object*/
  double **T = T_from_python(dataObj, n, k); /*Convert python T to C T*/
  spk_python_2(T, keys, n, k); /*Finish spk and kmeans algorithm with C*/
  return Py_BuildValue("");
}
/*Convert T from python to to double** T in C*/
double **T_from_python(PyObject *dataObj, int n, int k)
{
  double **T = (double **)malloc(sizeof(double *) * n);
  if (T == NULL)
  {
    exit_run();
  }
  for (i = 0; i < n; i++)
  {
    T[i] = (double *)malloc(sizeof(double) * k);
    if (T[i] == NULL)
    {
      exit_run();
    }
  }
  for (i = 0; i < n; i++)
  {
    PyObject *Pointobj = PyList_GetItem(dataObj, i);

    for (j = 0; j < k; j++)
    {
      PyObject *coor_point = PyList_GetItem(Pointobj, j);

      double coor_po = PyFloat_AS_DOUBLE(coor_point);
      T[i][j] = coor_po;
    }
  }

  return T;
}
/*First part of spk in C- create laplacian and perform Jacobi*/
static PyObject *spk_python_1(int k, char *filename)
{
  /*Regular C functions*/
  int *d_n = get_d_n(filename);
  n = d_n[1];
  if (k >= n)
  {
    printf("Invalid Input!");
    return Py_BuildValue("");
  }
  d = d_n[0];
  data = build_arr(d_n[0], d_n[1], filename); /* build datapoints
   from the input file*/
  W = build_w(data, n);
  D = build_D(W, n);
  free_point(data, n);
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
  T = getT(U, n, k);             /* Get the matrix T as the datapoints
   - K means will be performed on T*/
  free_mat(U, n);
  free_mat(V, n);
  free(eigengap);
  free(tup);
  return T_Python(n, T, k); /*T is sent back to python as PyObject*/
}

/*Perform the second part of spk, Kmeans with the selected clusters*/
void spk_python_2(double **T, double *keys, int n, int k)
{
  datapoints = build_DataPoints(T, n, k);
  clusters = create_py_clusters(datapoints, k, k, keys);
  free(keys);
  for (i = 0; i < n; i++)
  {
    free(T[i]);
  }
  free(T);
  k_means(datapoints, clusters, k, k, n, 300);
  for (i = 0; i < n; i++)
  {
    free(datapoints[i].coords);
  }
  free(datapoints);

  free(clusters);
}

/*Convert python object  keys to C keys*/
double *get_keys_from_py(PyObject *listObj, int k)
{
  double *keys;
  PyObject *key_py;
  keys = (double *)malloc(sizeof(double) * k);
  if (keys == NULL)
  {
    exit_run();
  }
  for (i = 0; i < k; i++)
  {
    key_py = PyList_GetItem(listObj, i);
    long key_long = PyLong_AsLong(key_py);
    int key = (int)key_long;
    keys[i] = key;
  }

  return keys;
}

/*Create clusters according to kmeans++ with the given keys*/
cluster *create_py_clusters(point *data_points, int k, int d, double *keys)
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
    int key = keys[i];
    clusters[i].cent = data_points[key];
    clusters[i].sum = sum_arr;
    clusters[i].size = 0.0;
  };

  return clusters;
}
/*T is sent back to python as PyObject*/
static PyObject *T_Python(int n, double **T, int k)
{
  PyObject *l_p;
  PyObject *l2;
  l_p = PyList_New(n);
  for (i = 0; i < n; i++)
  {
    /* list of coordinates of point in the clusters*/
    l2 = PyList_New(k);
    for (j = 0; j < k; j++)
    {
      PyList_SetItem(l2, j, PyFloat_FromDouble(T[i][j]));
    }
    free(T[i]);

    PyList_SetItem(l_p, i, l2);
  }
  free(T);

  return l_p;
}

/* configurations */
static PyMethodDef _kmeansMethods[] = {
    {"standard_interface",
     (PyCFunction)standard_interface,
     METH_VARARGS,
     PyDoc_STR("func for regular")},
    {"spk_a",
     (PyCFunction)spk_a,
     METH_VARARGS,
     PyDoc_STR("func 1 for spk to send T ")},
    {"spk_b",
     (PyCFunction)spk_b,
     METH_VARARGS,
     PyDoc_STR("func 2 for spk to finish")},
    {NULL, NULL, 0, NULL}

};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "spk",
    NULL,
    -1,
    _kmeansMethods};

PyMODINIT_FUNC
PyInit_spk(void)
{
  return PyModule_Create(&_moduledef);
}