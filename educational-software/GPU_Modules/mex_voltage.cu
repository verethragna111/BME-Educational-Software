#include "mex.h"
#include "cuda_runtime.h"

extern void voltage_clamp_kernel(
    float* Im, float* INa, float* IK,
    float V_clamp0, float V_clamp1,
    float E_NA, float E_K,
    int idel, int idel2, float DT, int len);

// MEX gateway
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 8)
        mexErrMsgTxt("Usage: mex_voltage(V0, V1, ENa, EK, idel, idel2, dt, len)");

    float V0 = (float)mxGetScalar(prhs[0]);
    float V1 = (float)mxGetScalar(prhs[1]);
    float ENa = (float)mxGetScalar(prhs[2]);
    float EK = (float)mxGetScalar(prhs[3]);
    int idel = (int)mxGetScalar(prhs[4]);
    int idel2 = (int)mxGetScalar(prhs[5]);
    float DT = (float)mxGetScalar(prhs[6]);
    int len = (int)mxGetScalar(prhs[7]);

    size_t bytes = len * sizeof(float);

    float *d_Im, *d_INa, *d_IK;
    cudaMalloc(&d_Im, bytes);
    cudaMalloc(&d_INa, bytes);
    cudaMalloc(&d_IK, bytes);

    dim3 threadsPerBlock(256);
    dim3 blocksPerGrid((len + 255) / 256);

    voltage_clamp_kernel<<<blocksPerGrid, threadsPerBlock>>>(
        d_Im, d_INa, d_IK, V0, V1, ENa, EK, idel, idel2, DT, len
    );

    // Output to MATLAB
    plhs[0] = mxCreateNumericMatrix(1, len, mxSINGLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, len, mxSINGLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, len, mxSINGLE_CLASS, mxREAL);

    cudaMemcpy(mxGetData(plhs[0]), d_Im, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(mxGetData(plhs[1]), d_INa, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(mxGetData(plhs[2]), d_IK, bytes, cudaMemcpyDeviceToHost);

    cudaFree(d_Im);
    cudaFree(d_INa);
    cudaFree(d_IK);
}
