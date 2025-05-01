#include "mex.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cmath> // For exp, pow, fabs

// Define device helper functions for alpha and beta rates
__device__ double alpha_m_d(double V) {
    if (fabs(25.0 - V) < 1e-9) { // Handle V close to 25.0 mV
        return 1.0; // Limit as V approaches 25.0
    }
    return 0.1 * (25.0 - V) / (exp((25.0 - V) / 10.0) - 1.0);
}

__device__ double beta_m_d(double V) {
    return 4.0 * exp(-V / 18.0);
}

__device__ double alpha_h_d(double V) {
    return 0.07 * exp(-V / 20.0);
}

__device__ double beta_h_d(double V) {
    return 1.0 / (exp((30.0 - V) / 10.0) + 1.0);
}

__device__ double alpha_n_d(double V) {
    if (fabs(10.0 - V) < 1e-9) { // Handle V close to 10.0 mV
        return 0.1; // Limit as V approaches 10.0
    }
    return 0.01 * (10.0 - V) / (exp((10.0 - V) / 10.0) - 1.0);
}

__device__ double beta_n_d(double V) {
    return 0.125 * exp(-V / 80.0);
}

// CUDA Kernel to perform the Hodgkin-Huxley calculations for each time step
__global__ void voltageBE_kernel(
    int num_time_steps,
    double DT,
    int tdel1_idx,
    int tdel2_idx,
    double V_REST,
    double G_NA,
    double G_K,
    double E_K,
    double E_NA,
    double V_clamp1,
    double V_clamp0,
    double m0,
    double h0,
    double n0,
    double* d_Im,
    double* d_INa,
    double* d_IK,
    double* d_g_Na,
    double* d_g_K)
{
    // Get the global thread index
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Ensure the thread index is within the bounds of the time steps
    if (i < num_time_steps) {

        double V;
        // Determine the voltage at the current time step based on indices
        if (i <= tdel1_idx) {
            V = V_REST;
        } else if (i >= tdel2_idx) {
            V = V_clamp1;
        } else {
            V = V_clamp0;
        }

        // Calculate steady-state values and time constants for the current voltage V
        double m_inf = alpha_m_d(V) / (alpha_m_d(V) + beta_m_d(V));
        double h_inf = alpha_h_d(V) / (alpha_h_d(V) + beta_h_d(V));
        double n_inf = alpha_n_d(V) / (alpha_n_d(V) + beta_n_d(V));

        double tau_m = 1.0 / (alpha_m_d(V) + beta_m_d(V));
        double tau_h = 1.0 / (alpha_h_d(V) + beta_h_d(V));
        double tau_n = 1.0 / (alpha_n_d(V) + beta_n_d(V));

        double M = m_inf - (m_inf - m0) * exp(-(i - tdel1_idx) * DT / tau_m);
        double H = h_inf - (h_inf - h0) * exp(-(i - tdel1_idx) * DT / tau_h);
        double N = n_inf - (n_inf - n0) * exp(-(i - tdel1_idx) * DT / tau_n);

        // Calculate conductances
        double gNa = G_NA * pow(M, 3) * H;
        double gK = G_K * pow(N, 4);

        // Calculate currents
        double Im = (gNa * (V - E_NA) + gK * (V - E_K));
        double INa = gNa * (V - E_NA);
        double IK = gK * (V - E_K);

        // Store results in global memory
        d_Im[i] = Im;
        d_INa[i] = INa;
        d_IK[i] = IK;
        d_g_Na[i] = gNa;
        d_g_K[i] = gK;
    }
}

// Define host helper functions for alpha and beta rates (for initial m, h, n calculation)
double alpha_m_h(double V) {
    if (fabs(25.0 - V) < 1e-9) {
        return 1.0;
    }
    return 0.1 * (25.0 - V) / (exp((25.0 - V) / 10.0) - 1.0);
}

double beta_m_h(double V) {
    return 4.0 * exp(-V / 18.0);
}

double alpha_h_h(double V) {
    return 0.07 * exp(-V / 20.0);
}

double beta_h_h(double V) {
    return 1.0 / (exp((30.0 - V) / 10.0) + 1.0);
}

double alpha_n_h(double V) {
    if (fabs(10.0 - V) < 1e-9) {
        return 0.1;
    }
    return 0.01 * (10.0 - V) / (exp((10.0 - V) / 10.0) - 1.0);
}

double beta_n_h(double V) {
    return 0.125 * exp(-V / 80.0);
}


// MEX function entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // --- Input Parsing ---
    // Expected inputs:
    // prhs[0]: V_clamp1 (double)
    // prhs[1]: Nae (double)
    // prhs[2]: tdel1 (double)
    // prhs[3]: tend (double)
    // prhs[4]: use_modified_flag (double, 0 or 1)
    // prhs[5]: V_clamp0 (double, required if use_modified_flag is 1)
    // prhs[6]: tdel2 (double, required if use_modified_flag is 1)

    if (nrhs < 5 || nrhs > 7) {
        mexErrMsgTxt("Incorrect number of input arguments.");
    }
    if (nlhs < 5 || nlhs > 6) {
         mexErrMsgTxt("Incorrect number of output arguments.");
    }


    double V_clamp1 = mxGetScalar(prhs[0]);
    double Nae = mxGetScalar(prhs[1]);
    double tdel1 = mxGetScalar(prhs[2]);
    double tend = mxGetScalar(prhs[3]);
    bool use_modified = (mxGetScalar(prhs[4]) != 0);

    double V_clamp0 = V_clamp1; // Default if not using modified
    double tdel2 = tend;       // Default if not using modified

    if (use_modified) {
        if (nrhs < 7) {
             mexErrMsgTxt("V_clamp0 and tdel2 must be provided if use_modified_flag is true.");
        }
        V_clamp0 = mxGetScalar(prhs[5]);
        tdel2 = mxGetScalar(prhs[6]);
    }

    // --- Constants and Initial Calculations ---
    const double DT = 0.05;
    const double V_REST = 0.0;
    const double G_NA = 120.0;	// mS/cm^2
    const double G_K = 25.0;	    // mS/cm^2
    // const double G_LEAK = 0.3;	// mS/cm^2 (Not used in the provided loop calculation)

    const double E_K = -8.0;         // mV
    // const double E_LEAK = 10.613;	// mV (Not used in the provided loop calculation)

    // Calculate E_NA on the host
    double E_NA = 25.0 * log(Nae / 45.0) + 60.0;

    // Calculate initial (steady-state) gating variables at V_REST on the host
    double m0 = alpha_m_h(V_REST) / (alpha_m_h(V_REST) + beta_m_h(V_REST));
    double h0 = alpha_h_h(V_REST) / (alpha_h_h(V_REST) + beta_h_h(V_REST));
    double n0 = alpha_n_h(V_REST) / (alpha_n_h(V_REST) + beta_n_h(V_REST));

    // Calculate time step indices
    int tdel1_idx = floor(tdel1 / DT);
    int tdel2_idx = floor(tdel2 / DT);
    int num_time_steps = floor(tend / DT) + 1;

    // Ensure indices are non-negative
    if (tdel1_idx < 0) tdel1_idx = 0;
    if (tdel2_idx < 0) tdel2_idx = 0;


    // --- CUDA Memory Allocation and Data Transfer ---
    double *d_Im, *d_INa, *d_IK, *d_g_Na, *d_g_K;
    size_t array_size = num_time_steps * sizeof(double);

    // Allocate device memory for output arrays
    cudaMalloc(&d_Im, array_size);
    cudaMalloc(&d_INa, array_size);
    cudaMalloc(&d_IK, array_size);
    cudaMalloc(&d_g_Na, array_size);
    cudaMalloc(&d_g_K, array_size);

    // Check for CUDA allocation errors
    if (cudaGetLastError() != cudaSuccess) {
        mexErrMsgTxt("CUDA memory allocation failed.");
    }

    // --- Kernel Configuration and Launch ---
    int block_size = 256; // Number of threads per block
    int grid_size = (num_time_steps + block_size - 1) / block_size; // Number of blocks

    // Launch the CUDA kernel
    voltageBE_kernel<<<grid_size, block_size>>>(
        num_time_steps,
        DT,
        tdel1_idx,
        tdel2_idx,
        V_REST,
        G_NA,
        G_K,
        E_K,
        E_NA,
        V_clamp1,
        V_clamp0,
        m0,
        h0,
        n0,
        d_Im,
        d_INa,
        d_IK,
        d_g_Na,
        d_g_K);

    // Check for CUDA kernel launch errors
    if (cudaGetLastError() != cudaSuccess) {
        mexErrMsgTxt("CUDA kernel launch failed.");
    }

    // Synchronize device to ensure kernel completes
    cudaDeviceSynchronize();
    if (cudaGetLastError() != cudaSuccess) {
        mexErrMsgTxt("CUDA device synchronization failed.");
    }

    // --- Data Transfer Back to Host ---
    double *h_Im, *h_INa, *h_IK, *h_g_Na, *h_g_K;

    // Create MATLAB output arrays
    plhs[0] = mxCreateDoubleMatrix(1, num_time_steps, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, num_time_steps, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, num_time_steps, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, num_time_steps, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, num_time_steps, mxREAL);

    // Get pointers to MATLAB output array data
    h_Im = mxGetPr(plhs[0]);
    h_INa = mxGetPr(plhs[1]);
    h_IK = mxGetPr(plhs[2]);
    h_g_Na = mxGetPr(plhs[3]);
    h_g_K = mxGetPr(plhs[4]);

    // Copy results from device to host (MATLAB arrays)
    cudaMemcpy(h_Im, d_Im, array_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_INa, d_INa, array_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_IK, d_IK, array_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_g_Na, d_g_Na, array_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_g_K, d_g_K, array_size, cudaMemcpyDeviceToHost);

    // Check for CUDA memcpy errors
     if (cudaGetLastError() != cudaSuccess) {
        mexErrMsgTxt("CUDA memcpy failed.");
    }


    // --- Output E_NA ---
    if (nlhs > 5) {
        plhs[5] = mxCreateDoubleScalar(E_NA);
    }

    // --- Cleanup ---
    // Free device memory
    cudaFree(d_Im);
    cudaFree(d_INa);
    cudaFree(d_IK);
    cudaFree(d_g_Na);
    cudaFree(d_g_K);

    // Check for CUDA free errors
     if (cudaGetLastError() != cudaSuccess) {
        mexErrMsgTxt("CUDA free failed.");
    }
}