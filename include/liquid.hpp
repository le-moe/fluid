
struct FluidCube {
    int size;
    float dt;
    float diff;
    float visc;
    
    float *s;
    float *density;
    
    float *Vx;
    float *Vy;
    // float *Vz;

    float *Vx0;
    float *Vy0;
    // float *Vz0;
};


FluidCube *FluidCubeCreate(int size, int diffusion, int viscosity, float dt);
void FluidCubeFree(FluidCube *cube);
void FluidCubeAddDensity(FluidCube *cube, int x, int y, float amount);
void FluidCubeAddVelocity(FluidCube *cube, int x, int y, float amountX, float amountY);
void FluidCubeStep(FluidCube *cube);
