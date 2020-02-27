// simulation constants
float MAX_DENSITY = 1000;
int N = 200;                 // how many boxes in each row and col - ther higher the better, but slower)
int ITER = 10;               // how many iteration to run Gauss-Seidel relaxation - the higher the more accurate, but slower :/
int RES = 3;                 // how big each box is in pixel - the smaller the better, but hard to see.
float DIFF = 0;            // diffuse rate, needs tuning
float VISC = 0.0000001f;              // viscocity rate, needs tuning
float dt = 0.1;          // timestep
Fluid fluid;
float t = 0;

// how much density to add to fluid at each mouse update
float delta = 500f; // if using trackpad, if using physical mouse, use 200f only

// pause the simulation
boolean paused = true;

// mouse prev position
float prevMouseX;
float prevMouseY;

// mouse velocity
PVector mouseVec;

void settings() {
  size((N-2)*RES, (N-2)*RES);
}

void setup() {
  fluid = new Fluid(N, DIFF, VISC, dt);
  prevMouseX = mouseX;
  prevMouseY = mouseY;
  mouseVec = new PVector();
}

void mouseDragged() {
  int x = int(map(mouseX, 0, width, 1, N-2));
  int y = int(map(mouseY, 0, height, 1, N-2));
  
  int startX = constrain(x-1, 1, N-2);
  int endX = constrain(x+1, 1, N-2);
  
  int startY = constrain(y-1, 1, N-1);
  int endY = constrain(y+1, 1, N-1);
  
  // add density based on mouse position
  for (int i = startX; i <= endX; i++) {
    for (int j = startY; j <= endY; j++) {
      fluid.addDensity(i, j, delta);
    }
  }
  
  // add velocity
  fluid.addVel(x, y, mouseVec.x, mouseVec.y);
}

void keyReleased() {
  if (key == 'v')
    paused = !paused;
}

void draw() {
  // setting the mouse vel
  mouseVec.x = (mouseX - prevMouseX) / 10;
  mouseVec.y = (mouseY - prevMouseY) / 10;
  prevMouseX = mouseX;
  prevMouseY = mouseY;
  
  background(0);
  if (!paused) {
    fluid.step(ITER);
    fluid.dissolve(-1);
  }
  fluid.render();
  surface.setTitle("FPS: " + round(frameRate));
}
