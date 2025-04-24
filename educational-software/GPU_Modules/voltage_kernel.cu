#include <cmath>

extern "C" 
__global__ void voltage_clamp_kernel(
    float* Im, float* INa, float* IK,
    float V_clamp0, float V_clamp1,
    float E_NA, float E_K,
    int idel, int idel2, float DT, int len)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= len) return;

    float V = 0.0f;
    if (i <= idel)
        V = 0.0f;
    else if (i >= idel2)
        V = V_clamp1;
    else
        V = V_clamp0;

    // Basic HH-like placeholders (simplified)
    float m = 0.05f;
    float h = 0.6f;
    float n = 0.32f;

    float gNa = 120.0f * powf(m, 3.0f) * h;
    float gK = 25.0f * powf(n, 4.0f);

    Im[i] = gNa * (V - E_NA) + gK * (V - E_K);
    INa[i] = gNa * (V - E_NA);
    IK[i] = gK * (V - E_K);
}
