class Fluid {
  int size;                    // size, or resolution of simulation
  float dt;                    // timestep for integration
  float diff;                  // diffuse rate (expand)
  float visc;                  // viscosity, fluid thickness

  float[] s;                   // density from previous timestep
  float[] density;             // density at current timestep

  float[] velX;                // vels at current timestep
  float[] velY;

  float[] velX0;               // vels at prev timestep
  float[] velY0;

  public Fluid(int size, float diff, float visc, float dt) {
    this.size = size;
    this.diff = diff;
    this.visc = visc;
    this.dt = dt;

    // initialize the arrays, we use 1D array
    // with a helper function to index into it for efficiency
    this.s = new float[size * size];
    this.density = new float[size * size];
    this.velX = new float[size * size];
    this.velY = new float[size * size];
    this.velX0 = new float[size * size];
    this.velY0 = new float[size * size];
  }

  public void addDensity(int row, int col, float amount) {
    if (amount >= 0 && density[IX(row, col)] < MAX_DENSITY)
      density[IX(row, col)] += amount;
    else if (density[IX(row, col)] > 0)
      density[IX(row, col)] += amount;
  }

  public void addVel(int row, int col, float amountX, float amountY) {
    velX[IX(row, col)] += amountX;
    velY[IX(row, col)] += amountY;
  }   

  public void diffuse(int b, float[] x, float[] x0, float diff, float dt, int iter) {
    float a = dt * diff * (size-2) * (size-2);
    linSolve(b, x, x0, a, 1+4*a, iter);
  } 

  public void project(float[] velocX, float[] velocY, float[] p, float[] div, int iter) {
    for (int j = 1; j < size - 1; j++) {
      for (int i = 1; i < size - 1; i++) {
        div[IX(i, j)] = -0.5f*(
          velocX[IX(i+1, j)]
          -velocX[IX(i-1, j)]
          +velocY[IX(i, j+1)]
          -velocY[IX(i, j-1)]
          )/size;
        p[IX(i, j)] = 0;
      }
    }
    setBound(0, div); 
    setBound(0, p);
    linSolve(0, p, div, 1, 4, iter);

    for (int j = 1; j < size - 1; j++) {
      for (int i = 1; i < size - 1; i++) {
        velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
          -p[IX(i-1, j)]) * size;
        velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
          -p[IX(i, j-1)]) * size;
      }
    }

    setBound(1, velocX);
    setBound(2, velocY);
  }

  public void advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float dt) {
    float i0, i1, j0, j1;

    float dtx = dt * (size - 2);
    float dty = dt * (size - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float Nfloat = size;
    float ifloat, jfloat;
    int i, j;

    for (j = 1, jfloat = 1; j < size - 1; j++, jfloat++) { 
      for (i = 1, ifloat = 1; i < size - 1; i++, ifloat++) {
        tmp1 = dtx * velocX[IX(i, j)];
        tmp2 = dty * velocY[IX(i, j)];
        x    = ifloat - tmp1; 
        y    = jfloat - tmp2;

        if (x < 0.5f) x = 0.5f; 
        if (x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
        i0 = floor(x); 
        i1 = i0 + 1.0f;
        if (y < 0.5f) y = 0.5f; 
        if (y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
        j0 = floor(y);
        j1 = j0 + 1.0f;

        s1 = x - i0; 
        s0 = 1.0f - s1; 
        t1 = y - j0; 
        t0 = 1.0f - t1;

        int i0i = int(i0);
        int i1i = int(i1);
        int j0i = int(j0);
        int j1i = int(j1);

        d[IX(i, j)] = 
          s0*(t0*d0[IX(i0i, j0i)]+t1*d0[IX(i0i, j1i)])+
          s1*(t0*d0[IX(i1i, j0i)]+t1*d0[IX(i1i, j1i)]);
      }
    }
    setBound(b, d);
  }

  public void setBound(int b, float[] x) {
    // set velocity for y
    for (int i = 1; i < size-1; i++) {
      x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
      x[IX(i, size-1)] = b == 2 ? -x[IX(i, size-2)] : x[IX(i, size-2)];
    }

    // set velocity for x
    for (int i = 1; i < size-1; i++) {
      x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
      x[IX(size-1, i)] = b == 1 ? -x[IX(size-2, i)] : x[IX(size-2, i)];
    }

    // set border vel and density
    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, size-1)] = 0.5 * (x[IX(1, size-1)] + x[IX(0, size-2)]);
    x[IX(size-1, 0)] = 0.5 * (x[IX(size-1, 1)] + x[IX(size-2, 0)]);
    x[IX(size-1, size-1)] = 0.5 * (x[IX(size-2, size-1)] + x[IX(size-1, size-2)]);
  }

  public void linSolve(int b, float[] x, float[] x0, float a, float c, int iter) {
    float cRecip = 1.0/c;
    for (int i = 0; i < iter; i++) {
      for (int k = 1; k < size-1; k++) {
        for (int j = 1; j < size-1; j++) {
          x[IX(j, k)] = (x0[IX(j, k)] + a * (x[IX(j+1, k)]
            + x[IX(j-1, k)]
            + x[IX(j, k+1)]
            + x[IX(j, k-1)])) * cRecip;
        }
      }
      setBound(b, x);
    }
  }

  public void step(int iter) {
    // this is the heart of the simulation
    diffuse(1, velX0, velX, visc, dt, iter);
    diffuse(2, velY0, velY, visc, dt, iter);

    project(velX0, velY0, velX, velY, iter);

    advect(1, velX, velX0, velX0, velY0, dt);
    advect(2, velY, velY0, velX0, velY0, dt);
    
    project(velX, velY, velX0, velY0, iter);
    
    diffuse(0, s, density, diff, dt, iter);
    advect(0, density, s, velX, velY, dt);
  }
  
  public void render() {
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        float col = map(density[IX(i, j)], 0, MAX_DENSITY, 0, 255.0);
        fill(col);
        noStroke();
        square((i-1)*RES, (j-1)*RES, RES);
      }
    }
  }
  
  public void dissolve(float rate) {
    assert rate < 0;
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        this.addDensity(i, j, rate);
      }
    }
  }

  public int IX(int x, int y) {
    x = constrain(x, 0, size-1);
    y = constrain(y, 0, size-1);
    return x + y*size;
  }
}
