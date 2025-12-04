typedef struct {
  float h1, h2, h3;
  int n1, n2, n3;
  float x_source, y_source, z_source;
  int niter, nfpi;
  float epsilon;
  float*** Vnmo;
  float*** V0;
  float*** eta;

  float*** T;
  float*** TT;
  float*** T0;
  float*** px0;
  float*** py0;
  float*** pz0;
  float*** tau;
  float*** px;
  float*** py;
  float*** pz;
  float*** rhs;

  int shotx[2], shoty[2], shotz[2];
} eikonal_t;
