void setup() { //<>//
  size(750, 750, P3D);
  surface.setTitle("Final Project - Smoke");
  //dens_randomize();
}

int N = 50; //smoke resolution
float kDiff = 0.01; //diffusion constant
float dens[][] = new float[N][N]; //density is from 0 - 100
float dens0[][] = new float[N][N];
float source[][] = new float[N][N];

void add_source(float dt) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      dens[i][j] += source[i][j] * dt;
    }
  }
}

void diffuse(float dt) {
  float a = dt * kDiff * N * N;
  
  for (int k = 0; k < 20; k++) { //Gauss-Sidel iterations
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        float sum = 0;
        if(i-1>0) sum += dens[i-1][j];
        if(i+1<N) sum += dens[i+1][j];
        if(j-1>0) sum += dens[i][j-1];
        if(j+1<N) sum += dens[i][j+1];
        dens[i][j] = (dens0[i][j] + a*(sum))/(1+4*a); 
      }
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

void update(float dt) {
  add_source(dt);
  //store prev dens
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      dens0[i][j] = dens[i][j];
    }
  }
  diffuse(dt);
  
}

boolean sourceMode = false; // if it's true, pressing the mouse will add a source in that position instead of just a little smoke
void mouseDragged() {
  if (mouseX < height && mouseX > 0 && mouseY > 0 && mouseY < height) { // make sure mouse is actually in the window first
    if(sourceMode) source[mouseX * N / height][mouseY * N / height] += 70;
    else dens[mouseX * N / height][mouseY * N / height] += 70;
  }
}

void mouseClicked() {
  if (mouseX < height && mouseX > 0 && mouseY > 0 && mouseY < height) { // make sure mouse is actually in the window first
    if(sourceMode) source[mouseX * N / height][mouseY * N / height] += 70;
    else dens[mouseX * N / height][mouseY * N / height] += 70;
  }
}

void keyPressed() {
  if(key == 'r') {
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
      float c = dens[i][j]/100 * 255;
      fill(c);
      rect((i * d), (j * d), d, d);
    }
  }
}
