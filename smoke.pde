void setup() { //<>//
  size(750, 750, P3D);
  surface.setTitle("Final Project - Smoke");
  //dens_randomize();

  //initialize velocity field
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      vx[i][j] = 0;
      vy[i][j] = -0.5;
    }
  }
}

int N = 50; //smoke resolution
float kDiff = 0.01; //diffusion constant
float dens[][] = new float[N][N]; //density is from 0 - 100
float dens0[][] = new float[N][N];
float source[][] = new float[N][N];
float vx[][] = new float[N][N];
float vy[][] = new float[N][N];

void add_source(float dt) {
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

void diffuse(float dt) {
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

void swap() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      float temp = dens[i][j];
      dens[i][j] = dens0[i][j];
      dens0[i][j] = temp;
    }
  }
}

void advect(float dt) {
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
      if (i1 < N && j1 < N) sum1 += t0*dens0[i1][j0]+t1*dens0[i1][j1];

      dens[i][j] = s0*sum0 + s1*sum1;
    }
  }
}

void update(float dt) {
  add_source(dt);
  swap();
  diffuse(dt);
  swap();
  advect(dt);
}

void mouseDens() {
  if (mouseX < height && mouseX > 0 && mouseY > 0 && mouseY < height) { // make sure mouse is actually in the window first
    if (sourceMode) source[mouseX * N / height][mouseY * N / height] += 70;
    else dens[mouseX * N / height][mouseY * N / height] += 70;
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

void draw() {

  //println("fps:", frameRate);
  background(0);
  noStroke();
  fill(255, 255, 255);

  float dt = 1/frameRate;
  update(dt);

  //draw grid contents
  rectMode(CORNER);
  float d = height/N;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      float c = (dens[i][j]/100) * 255;
      fill(c);
      rect((i * d), (j * d), d, d);
    }
  }
}
