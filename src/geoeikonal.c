#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <float.h>
#include <math.h>



// Min-heap A min-heap is used to track points at the edge of the front.  Only the index of each
// point is stored in the heap, the values are stored in the distance matrix.  The status array
// records whether a point is frozen so its value is known (status=-2), is far from the front and as
// yet unvisited (status=-1), or if it is a trial point on the front that may still be updated
// (status>=0). For trial points, the status is the index of the point in the heap, and is used to
// reheapify the heap when the point is updated.
typedef struct {
  int* point;
  int size;
  int capacity;
} MinHeap;

// Reheapify after the point at position k in the heap's valuehas decreased
static void reheapUp(MinHeap* heap, int k, const double *D, int* status, int point) {
  
  // Bubble up to maintain heap property
  const double value = D[point];
  while(k > 0) {
    int p = (k-1)/2;
    if(D[heap->point[p]] <= value) break;
    heap->point[k] = heap->point[p];
    if(status[heap->point[k]] >= 0)
      status[heap->point[k]] = k;
    k = p;
  }
  heap->point[k] = point;
  status[point] = k;
}

// Push a point onto the heap and mark it as a trial point
static void push(MinHeap* heap, const double* D, int* status, int point) {
  
  // Resize if at capacity
  if(heap->size == heap->capacity) {
    heap->capacity *= 2;
    heap->point = R_Realloc(heap->point, heap->capacity, int);
  }
  
  // Insert at end of array
  int k = heap->size++;
  heap->point[k] = point;

  // Reorder the heap
  reheapUp(heap, k, D, status, point);

}

// Reheapify after a points value has decreased
static void decrease(MinHeap* heap, const double *D, int* status, int point) {
  // Find point in the heap
  int k = status[point];

  // Reorder the heap
  reheapUp(heap, k, D, status, point);
}


// Pop the point with the smallest value from the heap
static int pop(MinHeap* heap, double *D, int* status) {
  int root = heap->point[0];
  int last = heap->point[--heap->size];

  // Bubble down to maintain heap property
  int p = 0, k = 0;
  while((p = 2*k+1) < heap->size) {
    if(p+1 < heap->size && D[heap->point[p+1]] < D[heap->point[p]]) p++;
    if(D[last] <= D[heap->point[p]]) break;
    heap->point[k] = heap->point[p];
    if(status[heap->point[k]] >= 0)
      status[heap->point[k]] = k;
    k = p;
  }
  heap->point[k] = last;
  status[last] = k;
  // Mark old root as frozen and return
  status[root] = -2;
  return root;  
}


static double solve_upwind(double vx, double hx, double vy, double hy, double f) {

  // If one of the points is infinite, return the other point plus the travel cost
  if (!R_finite(vx)) return vy+f*hy;
  if (!R_finite(vy)) return vx+f*hx;

  double hx2 = hx*hx;
  double hy2 = hy*hy;
  double a = 1.0/hx2 + 1.0/hy2;
  double b = -(vx/hx2+vy/hy2);
  double c = (vx*vx)/hx2+(vy*vy)/hy2-f*f;
  double d = b*b-a*c;
  
  if(d < 0.0) return fmin(vx+hx*f, vy+hy*f);
  double value = (sqrt(d)-b)/a;
  if(value < fmax(vx, vy)) return fmin(vx+hx*f, vy+hy*f);
  return value;
}


// Solve the eikonal equation on a regular grid using Sethian's upwinding fast marching
// method.
//
// nx: number of grid points in the x direction
// ny: number of grid points in the y direction
// Cost: A matrix of (positive) travel costs
// Dist: A matrix to store the solution
// hx: a vector of the grid spacing in the x direction
// hy: the grid spacing in the y direction
//
// The Cost matrix should be initliaised to 0 at source point, NA or negative at barrier points,
// and a postive cost at points where the solution is desired.
static void C_fastmarch(int nx, int ny, const double *Cost, double* Dist, 
                        const double *hx, const double hy) {
  // Total number of grid cells
  const int ntot = nx*ny;

  // Initialize min heap
  MinHeap heap;
  heap.size = 0;
  heap.capacity = (nx>ny)?nx:ny;
  heap.point = R_Calloc(heap.capacity, int);

  // Status array
  int* status = (int*)R_alloc(ntot, sizeof(int));

  // Initialize Dist, status and add sources to heap
  for(int k=0; k<ntot; k++) {
    // Set as far (unvisited)
    status[k] = -1;
    Dist[k] = R_PosInf;
    // Source point - add to heap and mark frozen
    if(Cost[k] == 0) {
      Dist[k] = 0;
      push(&heap,Dist,status,k);
      status[k] = -2;
      continue;
    }
    // Barrier point - mark frozen
    if(!R_finite(Cost[k]) || Cost[k] < 0) {
      Dist[k] = NA_REAL;
      status[k] = -2; 
    } 
  }

  // Neighbour offsets (left, right, up, down)
  const int di[4] = {-1, 1, 0, 0};
  const int dj[4] = {0, 0, -1, 1};

  // Main fast marching loop
  while(heap.size > 0) {
    // Extract point with smallest current value
    const int k = pop(&heap,Dist,status);
    const int i = k%nx;
    const int j = k/nx;
    
    // Scan the four neighbours
    for(int n=0; n<4; n++) {
      const int ni = i+di[n];
      const int nj = j+dj[n];
      const int nk = ni+nj*nx;
      
      // Skip if point is outside grid
      if(ni < 0 || ni >= nx || nj < 0 || nj >= ny) continue;
      
      // Skip if point frozen or barrier
      if(status[nk]==-2) continue;
      
      // Get minimum value of x neighbours
      double vx = R_PosInf;
      if(ni > 0 && R_finite(Dist[nk-1])) vx = Dist[nk-1];
      if(ni < nx-1 && R_finite(Dist[nk+1])) vx = fmin(vx, Dist[nk+1]);

      // Get minimum value of y neighbours
      double vy = R_PosInf;
      if(nj > 0 && R_finite(Dist[nk-nx])) vy = Dist[nk-nx];
      if(nj < ny-1 && R_finite(Dist[nk+nx])) vy = fmin(vy, Dist[nk+nx]);
      
      // Solve upwind quadratic
      double new_value = solve_upwind(vx, hx[nj], vy, hy, Cost[nk]);
      
      // Update if point is far or new value is smaller
      if(status[nk]==-1 || new_value < Dist[nk]) {
        Dist[nk] = new_value;
        // If point is already in heap, reheapify, else push
        if(status[nk]>=0)
          decrease(&heap,Dist,status,nk);
        else
          push(&heap,Dist,status,nk);
      }
    }
  }
  
  // Free min heap
  R_Free(heap.point);
}


// Solve the eikonal equation on a regular grid using Sethian's upwinding fast marching
// method.
//
// nx: number of grid points in the x direction
// ny: number of grid points in the y direction
// Cost: A matrix of (positive) travel costs
// Dist: A matrix to store the solution
// hx: a vector of the grid spacing in the x direction
// hy: the grid spacing in the y direction
// wrapx: should the grid wrap in the x direction
//
// The Cost matrix should be initliaised to 0 at source point, NA or negative at barrier points,
// and a postive cost at points where the solution is desired.
static void C_fastmarchp(int nx, int ny, const double *Cost, double* Dist, 
                         const double *hx, const double hy, const int wrapx) {
  // Total number of grid cells
  const int ntot = nx*ny;

  // Initialize min heap
  MinHeap heap;
  heap.size = 0;
  heap.capacity = (nx>ny)?nx:ny;
  heap.point = R_Calloc(heap.capacity, int);

  // Status array
  int* status = (int*)R_alloc(ntot, sizeof(int));

  // Initialize Dist, status and add sources to heap
  for(int k=0; k<ntot; k++) {
    // Set as far (unvisited)
    status[k] = -1;
    Dist[k] = R_PosInf;
    // Source point - add to heap and mark frozen
    if(Cost[k] == 0) {
      Dist[k] = 0;
      push(&heap,Dist,status,k);
      status[k] = -2;
      continue;
    }
    // Barrier point - mark frozen
    if(!R_finite(Cost[k]) || Cost[k] < 0) {
      Dist[k] = NA_REAL;
      status[k] = -2; 
    }
  }

  // Main fast marching loop
  while(heap.size > 0) {
    // Extract point with smallest current value
    const int k = pop(&heap,Dist,status);
    const int i = k%nx;
    const int j = k/nx;

    // Neighbour offsets (left, right, up, down)
    int di[4] = {-1, 1, 0, 0};
    int dj[4] = {0, 0, -1, 1};
    if(wrapx) {
      if(i==0) di[0] = nx-1;
      if(i==nx-1) di[1] = 1-nx;
    }
 
    // Scan the four neighbours
    for(int n=0; n<4; n++) {
      const int ni = i+di[n];
      const int nj = j+dj[n];
      const int nk = ni+nj*nx;
      
      // Skip if point is outside grid
      if(ni < 0 || ni >= nx || nj < 0 || nj >= ny) continue;
      
      // Skip if point frozen or barrier
      if(status[nk]==-2) continue;
      
      // Get minimum value of x neighbours
      double vx = R_PosInf;
      if(wrapx) {
        int offset = (ni>0) ? -1 : (nx-1);
        if(R_finite(Dist[nk+offset])) vx = Dist[nk+offset];
        offset = (ni<nx-1) ? 1 : (1-nx);
        if(R_finite(Dist[nk+offset])) vx = fmin(vx, Dist[nk+offset]);
      } else {
        if(ni > 0 && R_finite(Dist[nk-1])) vx = Dist[nk-1];
        if(ni < nx-1 && R_finite(Dist[nk+1])) vx = fmin(vx, Dist[nk+1]);
      }

      // Get minimum value of y neighbours
      double vy = R_PosInf;
      if(nj > 0 && R_finite(Dist[nk-nx])) vy = Dist[nk-nx];
      if(nj < ny-1 && R_finite(Dist[nk+nx])) vy = fmin(vy, Dist[nk+nx]);
      
      // Solve upwind quadratic
      double new_value = solve_upwind(vx, hx[nj], vy, hy, Cost[nk]);
      
      // Update if point is far or new value is smaller
      if(status[nk]==-1 || new_value < Dist[nk]) {
        Dist[nk] = new_value;
        // If point is already in heap, reheapify, else push
        if(status[nk]>=0)
          decrease(&heap,Dist,status,nk);
        else
          push(&heap,Dist,status,nk);
      }
    }
  }
  
  // Free min heap
  R_Free(heap.point);
}


// Calculate the upwind gradient of a solution to the eikonal equation
//
// nx: number of grid points in the x direction
// ny: number of grid points in the y direction
// Dist: A solution of the eikonal equation
// Gx: A matrix to store the gradient in the x direction
// Gy: A matrix to store the gradient in the y direction
// hx: a vector of the grid spacing in the x direction
// hy: the grid spacing in the y direction
//
// The Dist matrix should be obtained from a call to C_fastmarch or C_fastmarchp
static void C_gradient(int nx, int ny, const double* Dist, double* Gx, double* Gy, 
                        const double *hx, const double hy) {
  // Total number of grid cells
  const int ntot = nx*ny;

  // Loop over grid
  for(int k=0; k<ntot; k++) {
    // Extract point with smallest current value
    const int i = k%nx;
    const int j = k/nx;
  
    double gx = NA_REAL;
    double gy = NA_REAL;
    // Skip if point is unreachable
    if(R_finite(Dist[k])) {

      // Upwind gradient in x direction
      int kl = k-1;
      int kr = k+1;
      int left = (i > 0 && R_finite(Dist[kl]));
      int right = (i < nx-1 && R_finite(Dist[kr]));
      if(left && right) {
        // Both neighbours are finite - use the smaller
        if(Dist[kl] < Dist[kr])
          gx = (Dist[k] < Dist[kl]) ? 0.0 : (Dist[k] - Dist[kl])/hx[j];
        else
          gx = (Dist[k] < Dist[kr]) ? 0.0 : (Dist[kr] - Dist[k])/hx[j];
      } else {
        // At most one neighbour is finite
        if(left) gx = (Dist[k] < Dist[kl]) ? 0.0 :(Dist[k] - Dist[kl])/hx[j];
        if(right) gx = (Dist[k] < Dist[kr]) ? 0.0 : (Dist[kr] - Dist[k])/hx[j];
      }

      // Upwind gradient in y direction
      int ku = k-nx;
      int kd = k+nx;
      int up = (j > 0 && R_finite(Dist[ku]));
      int down = (j < ny-1 && R_finite(Dist[kd]));
      if(up && down) {
        // Both neighbours are finite - use the smaller
        if(Dist[ku] < Dist[kd])
          gy = (Dist[k] < Dist[ku]) ? 0.0 : (Dist[k] - Dist[ku])/hy;
        else
          gy = (Dist[k] < Dist[kd]) ? 0.0 : (Dist[kd] - Dist[k])/hy;
      } else {
        // At most one neighbour is finite
        if(up) gy = (Dist[k] < Dist[ku]) ? 0.0 : (Dist[k] - Dist[ku])/hy;
        if(down) gy = (Dist[k] < Dist[kd]) ? 0.0 : (Dist[kd] - Dist[k])/hy;
      }
    }
    Gx[k] = gx;
    Gy[k] = gy;
  }
}




// Calculate the upwind gradient of a solution to the eikonal equation
//
// nx: number of grid points in the x direction
// ny: number of grid points in the y direction
// Dist: A solution of the eikonal equation
// Gx: A matrix to store the gradient in the x direction
// Gy: A matrix to store the gradient in the y direction
// hx: a vector of the grid spacing in the x direction
// hy: the grid spacing in the y direction
// wrapx: should the grid wrap in the x direction
//
// The Dist matrix should be obtained from a call to C_fastmarch or C_fastmarchp
static void C_gradientp(int nx, int ny, const double* Dist, double* Gx, double* Gy, 
                        const double *hx, const double hy, const int wrapx) {
  // Total number of grid cells
  const int ntot = nx*ny;

  // Loop over grid
  for(int k=0; k<ntot; k++) {
    // Extract point with smallest current value
    const int i = k%nx;
    const int j = k/nx;
  
    double gx = NA_REAL;
    double gy = NA_REAL;
    // Skip if point is unreachable
    if(R_finite(Dist[k])) {

      // Upwind gradient in x direction
      int kl = k-1;
      int kr = k+1;
      int left = 0, right = 0;
      if(wrapx) {
        if(i==0) kl += nx;
        if(i==nx-1) kr -= nx;
        left = R_finite(Dist[kl]);
        right = R_finite(Dist[kr]);
      } else {
        left = (i > 0 && R_finite(Dist[kl]));
        right = (i < nx-1 && R_finite(Dist[kr]));
      }
      if(left && right) {
        // Both neighbours are finite - use the smaller
        if(Dist[kl] < Dist[kr])
          gx = (Dist[k] < Dist[kl]) ? 0.0 : (Dist[k] - Dist[kl])/hx[j];
        else
          gx = (Dist[k] < Dist[kr]) ? 0.0 : (Dist[kr] - Dist[k])/hx[j];
      } else {
        // At most one neighbour is finite
        if(left) gx = (Dist[k] < Dist[kl]) ? 0.0 : (Dist[k] - Dist[kl])/hx[j];
        if(right) gx = (Dist[k] < Dist[kr]) ? 0.0 : (Dist[kr] - Dist[k])/hx[j];
      }

      // Upwind gradient in y direction
      int ku = k-nx;
      int kd = k+nx;
      int up = (j > 0 && R_finite(Dist[ku]));
      int down = (j < ny-1 && R_finite(Dist[kd]));
      if(up && down) {
        // Both neighbours are finite - use the smaller
        if(Dist[ku] < Dist[kd])
          gy = (Dist[k] < Dist[ku]) ? 0.0 : (Dist[k] - Dist[ku])/hy;
        else
          gy = (Dist[k] < Dist[kd]) ? 0.0 : (Dist[kd] - Dist[k])/hy;
      } else {
        // At most one neighbour is finite
        if(up) gy = (Dist[k] < Dist[ku]) ? 0.0 : (Dist[k] - Dist[ku])/hy;
        if(down) gy = (Dist[k] < Dist[kd]) ? 0.0 : (Dist[kd] - Dist[k])/hy;
      }
    }
    Gx[k] = gx;
    Gy[k] = gy;
  }
}


// Solve the eikonal equation on a longitude-latitude grid using Sethian's upwinding fast marching
// method.
//
// cost_matrix: A matrix of (positive) travel costs
// extent_vec: geographic extent of the grid c(lon_min, lon_max, lat_min, lat_max)
// result_matrix: A matrix to store the solution
//
// The Cost matrix should be initliaised to 0 at source point, NA or negative at barrier points,
// and a postive cost at points where the solution is desired.
//
// Returns the solution matrix
SEXP C_geoeikonal(SEXP cost_matrix, SEXP extent_vec, SEXP result_matrix) {

  // Get dimensions and data
  const int nrow = Rf_nrows(cost_matrix);
  const int ncol = Rf_ncols(cost_matrix);
  const double* Cost = REAL(cost_matrix);
  const double* extent = REAL(extent_vec);
  double* Dist = REAL(result_matrix);

  // The (constant) column spacing
  const double R = 6378137;
  double dlat = M_PI/180.0*(extent[3]-extent[2])/ncol;
  double hlat = R*dlat;

  // The row spacing in each column
  double* hlon = (double*)R_alloc(ncol, sizeof(double));
  double dlon = M_PI/180.0*(extent[1]-extent[0])/nrow;
  for(int j=0; j<ncol; j++) {
    double lat = M_PI/180.0*extent[3] - (j+0.5)*dlat;
    //dlon[j] = R*acos(sin(lat)*sin(lat)+cos(lat)*cos(lat)*cos(dlon));
    hlon[j] = R*cos(lat)*dlon;
  }
  int wrapx = (extent[1] - extent[0] == 360);

  // Solve the eikonal equation
  C_fastmarchp(nrow, ncol, Cost, Dist, hlon, hlat, wrapx);

  return result_matrix;
}


// Solve the eikonal equation on a regular grid using Sethian's upwinding fast marching
// method.
//
// cost_matrix: A matrix of (positive) travel costs
// hrow: vector of row spacings in each column
// hcol: vector of column spacings in each row
// result_matrix: A matrix to store the solution
//
// The Cost matrix should be initliaised to 0 at source point, NA or negative at barrier points,
// and a postive cost at points where the solution is desired.
//
// Returns the solution matrix
SEXP C_grdeikonal(SEXP cost_matrix, SEXP hrow_vec, SEXP hcol_vec, SEXP result_matrix) {

  // Get dimensions and data
  const int nrow = Rf_nrows(cost_matrix);
  const int ncol = Rf_ncols(cost_matrix);
  const double* Cost = REAL(cost_matrix);
  const double* hrow0 = REAL(hrow_vec);
  const double* hcol0 = REAL(hcol_vec);
  double* Dist = REAL(result_matrix);


  // The (constant) column spacing
  double hcol = hcol0[0];

  // The row spacing in each column
  double* hrow = (double*)R_alloc(ncol, sizeof(double));
  for(int j=0; j<ncol; j++) hrow[j] = hrow0[0];

  // Solve the eikonal equation
  C_fastmarch(nrow, ncol, Cost, Dist, hrow, hcol);

  return result_matrix;
}


// Calculate the gradient of a solution to the eikonal equation on a longitude-latitude grid
//
// dist_matrix: A matrix of distances
// extent_vec: geographic extent of the grid c(lon_min, lon_max, lat_min, lat_max)
// gradx_matrix: A matrix to store the gradient in the x direction
// grady_matrix: A matrix to store the gradient in the y direction
//
// The Dist matrix should be obtained from a call to C_geoeikonal.
//
// Returns the components of the gradient as a two element list
SEXP C_geogradient(SEXP dist_matrix, SEXP extent_vec, SEXP gradx_matrix, SEXP grady_matrix) {

  // Get dimensions and data
  const int nrow = Rf_nrows(dist_matrix);
  const int ncol = Rf_ncols(dist_matrix);
  const double* extent = REAL(extent_vec);
  const double* Dist = REAL(dist_matrix);
  double* Gx = REAL(gradx_matrix);
  double* Gy = REAL(grady_matrix);

  // The (constant) column spacing
  const double R = 6378137;
  double dlat = M_PI/180.0*(extent[3]-extent[2])/ncol;
  double hlat = R*dlat;

  // The row spacing in each column
  double* hlon = (double*)R_alloc(ncol, sizeof(double));
  double dlon = M_PI/180.0*(extent[1]-extent[0])/nrow;
  for(int j=0; j<ncol; j++) {
    double lat = M_PI/180.0*extent[3] - (j+0.5)*dlat;
    //dlon[j] = R*acos(sin(lat)*sin(lat)+cos(lat)*cos(lat)*cos(dlon));
    hlon[j] = R*cos(lat)*dlon;
  }
  int wrapx = (extent[1] - extent[0] == 360);

  // Calculate the gradient
  C_gradientp(nrow, ncol, Dist, Gx, Gy, hlon, hlat, wrapx);

  // Return components as a list
  SEXP result_list = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result_list, 0, gradx_matrix);
  SET_VECTOR_ELT(result_list, 1, grady_matrix);

  // Set names for the list
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("lon"));
  SET_STRING_ELT(names, 1, Rf_mkChar("lat"));
  Rf_setAttrib(result_list, R_NamesSymbol, names);

  UNPROTECT(2);
  return result_list;
}


// Calculate the gradient of a solution to the eikonal equation on a regular grid
//
// dist_matrix: A matrix of distances
// hrow_vec: vector of row spacings in each column
// hcol_vec: vector of column spacings in each row
// gradx_matrix: A matrix to store the gradient in the x direction
// grady_matrix: A matrix to store the gradient in the y direction
//
// The Dist matrix should be obtained from a call to C_grdeikonal.
//
// Returns the components of the gradient as a two element list
SEXP C_grdgradient(SEXP dist_matrix, SEXP hrow_vec, SEXP hcol_vec, SEXP gradx_matrix, SEXP grady_matrix) {

  // Get dimensions and data
  const int nrow = Rf_nrows(dist_matrix);
  const int ncol = Rf_ncols(dist_matrix);
  const double* hrow0 = REAL(hrow_vec);
  const double* hcol0 = REAL(hcol_vec);
  const double* Dist = REAL(dist_matrix);
  double* Gx = REAL(gradx_matrix);
  double* Gy = REAL(grady_matrix);

  // The (constant) column spacing
  double hcol = hcol0[0];

  // The row spacing in each column
  double* hrow = (double*)R_alloc(ncol, sizeof(double));
  for(int j=0; j<ncol; j++) hrow[j] = hrow0[0];

  // Calculate the gradient
  C_gradient(nrow, ncol, Dist, Gx, Gy, hrow, hcol);

  // Return components as a list
  SEXP result_list = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result_list, 0, gradx_matrix);
  SET_VECTOR_ELT(result_list, 1, grady_matrix);

  // Set names for the list
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("x"));
  SET_STRING_ELT(names, 1, Rf_mkChar("y"));
  Rf_setAttrib(result_list, R_NamesSymbol, names);

  UNPROTECT(2);
  return result_list;
}



static const R_CallMethodDef CallEntries[] = {
    {"C_geoeikonal", (DL_FUNC) &C_geoeikonal, 3},
    {"C_grdeikonal", (DL_FUNC) &C_grdeikonal, 4},
    {"C_geogradient", (DL_FUNC) &C_geogradient, 4},
    {"C_grdgradient", (DL_FUNC) &C_grdgradient, 5},
    {NULL, NULL, 0}
};

void R_init_geoeikonal(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
} 