void setup() { //<>//
  size(750, 750, P3D);
  surface.setTitle("Final Project - Smoke");
  dens_randomize();

  init_vel();
  
}

int N = 100; //smoke resolution
float kDiff = 0.001; //diffusion constant
float visc = 0; //viscosity
float dens[][] = new float[N][N]; //density is from 0 - 100
float dens0[][] = new float[N][N];
float source[][] = new float[N][N];
float vx[][] = new float[N][N];
float vy[][] = new float[N][N];
float vx0[][] = new float[N][N];
float vy0[][] = new float[N][N];

//TODO - dens and source should be clamped to some max val
void add_source(float dens[][], float source[][], float dt) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      dens[i][j] += source[i][j] * dt;
    }
  }
}

void dens_randomize() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      dens[i][j] = random(100);
    }
  }
}

void init_vel() {
  //initialize velocity field
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      vx[i][j] = -random(-0.5, 0.5);
      vy[i][j] = -0.5;
    }
  }
}

void diffuse(float dens[][], float dens0[][], float kDiff, float dt) {
  float a = dt * kDiff * N * N;

  for (int k = 0; k < 20; k++) { //Gauss-Sidel iterations
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        float sum = 0;
        if (i-1>0) sum += dens[i-1][j];
        if (i+1<N) sum += dens[i+1][j];
        if (j-1>0) sum += dens[i][j-1];
        if (j+1<N) sum += dens[i][j+1];
        dens[i][j] = (dens0[i][j] + a*(sum))/(1+4*a);
      }
    }
  }
}

void swap(float x0[][], float x[][]) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      float temp = x[i][j];
      x[i][j] = x0[i][j];
      x0[i][j] = temp;
    }
  }
}

void advect(float dens[][], float dens0[][], float vx[][], float vy[][], float dt) {
  float dt0 = dt*N;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      float x = i-dt0*vx[i][j];
      float y = j-dt0*vy[i][j];

      if (x<0.5) x=0.5;
      if (x>N+0.5) x=N+0.5;

      int i0=int(x);
      int i1=i0+1;

      if (y<0.5) y=0.5;
      if (y>N+0.5) y=N+0.5;

      int j0=int(y);
      int j1=j0+1;

      float s1 = x-i0;
      float s0 = 1-s1;
      float t1 = y-j0;
      float t0 = 1-t1;

      float sum0 = 0;
      if (i0 < N && j0 < N) sum0 += t0*dens0[i0][j0];
      if (i0 < N && j1 < N) sum0 += t1*dens0[i0][j1];

      float sum1 = 0;
      if (i1 < N && j0 < N) sum1 += t0*dens0[i1][j0];
      if (i1 < N && j1 < N) sum1 += t1*dens0[i1][j1];

      dens[i][j] = s0*sum0 + s1*sum1;
    }
  }
}

//TODO:
//void set_bnd (int b, float x[][]) {
//  int i;
//  for ( i=1 ; i<N ; i++ ) {
//  x[IX(0 ,i)] = b==1 ? –x[IX(1,i)] : x[IX(1,i)];
//  x[IX(N+1,i)] = b==1 ? –x[IX(N,i)] : x[IX(N,i)];
//  x[IX(i,0 )] = b==2 ? –x[IX(i,1)] : x[IX(i,1)];
//  x[IX(i,N+1)] = b==2 ? –x[IX(i,N)] : x[IX(i,N)];
//  }
//  x[IX(0 ,0 )] = 0.5*(x[IX(1,0 )]+x[IX(0 ,1)]);
//  x[IX(0 ,N+1)] = 0.5*(x[IX(1,N+1)]+x[IX(0 ,N )]);
//  x[IX(N+1,0 )] = 0.5*(x[IX(N,0 )]+x[IX(N+1,1)]);
//  x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N )]);
//}

void project(float u[][], float v[][], float p[][], float div[][]) {
  float h = 1.0/N;
  float sum = 0;

  for (int i = 1; i < N; i++) {
    for (int j = 1; j < N; j++) {
      sum = 0;
      if (i-1>0) sum -= u[i-1][j];
      if (i+1<N) sum += u[i+1][j];
      if (j-1>0) sum -= v[i][j-1];
      if (j+1<N) sum += v[i][j+1];
      div[i][j] = -0.5*h*(sum);
      p[i][j] = 0;
    }
  }
  //set_bnd(0,div); set_bnd(0,p)

  for (int k = 0; k < 20; k++) { //Gauss-Sidel iterations
    for (int i = 1; i < N; i++) {
      for (int j = 1; j < N; j++) {
        sum = 0;
        if (i-1>0) sum += p[i-1][j];
        if (i+1<N) sum += p[i+1][j];
        if (j-1>0) sum += p[i][j-1];
        if (j+1<N) sum += p[i][j+1];
        p[i][j] = (div[i][j] + sum)/4;
      }
    }
  }
  //set_bnd(0,p)

  for (int i = 1; i < N; i++) {
    for (int j = 1; j < N; j++) {
      sum = 0;
      if (i-1>0) sum -= p[i-1][j];
      if (i+1<N) sum += p[i+1][j];
      u[i][j] -= 0.5*(sum)/h;

      sum = 0;
      if (j-1>0) sum -= p[i][j-1];
      if (j+1<N) sum += p[i][j+1];
      v[i][j] -= 0.5*(sum)/h;
    }
  }
  //set_bnd(1,u);set_bnd(2,v)
}

void update_dens(float dt) {
  add_source(dens, source, dt);
  swap(dens0, dens);
  diffuse(dens, dens0, kDiff, dt);
  swap(dens0, dens);
  advect(dens, dens0, vx, vy, dt);
}

void update_vel(float dt) {
  add_source(vx, vx0, dt);
  add_source(vy, vy0, dt);
  
  swap(vx0, vx);
  diffuse(vx, vx0, visc, dt);
  
  swap(vy0, vy);
  diffuse(vy, vy0, visc, dt);
  
  project(vx, vy, vx0, vy0);
  
  swap(vx0, vx);
  swap(vy0, vy);
  
  advect(vx, vx0, vx0, vy0, dt);
  advect(vy, vy0, vx0, vy0, dt);
  
  project(vx, vy, vx0, vy0);
}

void mouseDens() {
  if (mouseX < height && mouseX > 0 && mouseY > 0 && mouseY < height) { // make sure mouse is actually in the window first
    int i = (int)mouseX * N / height;
    int j = (int)mouseY * N / height;
    if (sourceMode) {
      source[i][j] += 70;
      if (i-1>0) source[i-1][j] += 70;
      if (i+1<N) source[i+1][j] += 70;
      if (j-1>0) source[i][j-1] += 70;
      if (j+1<N) source[i][j+1] += 70;
    } else {
      dens[i][j] += 70;
      if (i-1>0) dens[i-1][j] += 70;
      if (i+1<N) dens[i+1][j] += 70;
      if (j-1>0) dens[i][j-1] += 70;
      if (j+1<N) dens[i][j+1] += 70;
    }
  }
}

boolean sourceMode = false; // if it's true, pressing the mouse will add a source in that position instead of just a little smoke
void mouseDragged() {
  mouseDens();
}

void mousePressed() {
  mouseDens();
}

void keyPressed() {
  if (key == 'r') {
    sourceMode = !sourceMode;
    println("source mode =", sourceMode);
  }
}

int counter = 0;
void draw() {

  println("fps:", frameRate);
  background(0);
  noStroke();
  fill(255, 255, 255);

  float dt = 1/frameRate;
  if(counter % 5 == 0) init_vel(); //TODO - this is a placeholder for now, not sure what get_from_UI does in the paper
  update_vel(dt);
  update_dens(dt);
  counter++;
  //draw grid contents
  //TODO - still bugged when N>50
  rectMode(CORNER);
  float d = height/N;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      float c = (dens[i][j]/100) * 255;
      if (c != 0) {
        fill(c);
        rect((i * d), (j * d), d, d);
      }

      c = (source[i][j]/500) * 255;
      if (c != 0) {
        fill(c, 50, 50, 255);
        rect((i * d), (j * d), d, d);
      }
    }
  }
}
