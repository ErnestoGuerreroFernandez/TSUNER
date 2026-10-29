// TSUNER Copyright (C) 2025 Ernesto Guerrero Fernandez
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "GPU_kernels.cuh"
#include "SedFlux.cuh"

__global__ void GetDeltaT(double* d_V, double* d_deltaTV, SimParams params)
{
    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;

    if ((i < params.npx) && (j < params.npy)) {
        int pos  = j * params.npx + i;
        int nvol = params.nvol;

        double h   = d_V[pos];
        double2 q;
        q.x = d_V[1 * nvol + pos];
        q.y = d_V[2 * nvol + pos];

        double2 u;
        u_calculator(h, q, params.epsilon_h, u);

        double cg = sqrt(G * h);
        double lx = fabs(u.x) + cg;
        double ly = fabs(u.y) + cg;
        double dt = params.CFL * 0.5 * min(params.dx / lx, params.dy / ly);

        d_deltaTV[pos] = dt;
    }
}

static __device__ inline __attribute__((unused)) void interpolate_timeSPL(double *CC, double *tCC, long unsigned int timeCC, double current_time, double &h1, double &zb1, double2 &q1, int npx, int npy, int pos) {

    double dt = tCC[1] - tCC[0];
    int t = int((current_time - tCC[0]) / dt);
    t *= (t > 0);
    if (t > timeCC - 1)
        t = timeCC - 1;
    int tp = t + 1;
    int posCC  = t * (2 * npx + 2 * npy) * 4 + pos * 4;
    int posCCp = tp * (2 * npx + 2 * npy) * 4 + pos * 4;

    double m, b;
    m = (CC[posCCp] - CC[posCC]) / (tCC[tp] - tCC[t]);
    b = CC[posCC] - m * tCC[t];
    h1 = m * current_time + b;
    h1 *= UNITS;

    m = (CC[posCCp + 3] - CC[posCC + 3]) / (tCC[tp] - tCC[t]);
    b = CC[posCC + 3] - m * tCC[t];
    zb1 = m * current_time + b;

    h1 += zb1;

    m = (CC[posCCp + 1] - CC[posCC + 1]) / (tCC[tp] - tCC[t]);
    b = CC[posCC + 1] - m * tCC[t];
    q1.x = m * current_time + b;
    q1.x *= UNITS * h1;

    m = (CC[posCCp + 2] - CC[posCC + 2]) / (tCC[tp] - tCC[t]);
    b = CC[posCC + 2] - m * tCC[t];
    q1.y = m * current_time + b;
    q1.y *= UNITS * h1;
}

__global__ void hmax_calculator(double* d_hmax, double *d_V, double epsilon_h, int npx, int npy) {
    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;
    if ((i < npx) && (j < npy)) {
        int pos = j * npx + i;
        if (d_V[pos] > epsilon_h) {
            d_hmax[pos] = max(d_hmax[pos], d_V[pos]);
        }
    }
}

__global__ void bmin_calculator(double* d_bmin, double *d_V, int npx, int npy) {
    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;
    if ((i < npx) && (j < npy)) {
        int pos  = j * npx + i;
        int nvol = npx * npy;
        d_bmin[pos] = max(d_bmin[pos], d_V[3 * nvol + pos]);
    }
}

__global__ void eta_max_calculator(double* d_eta_max, double *d_V, double epsilon_h, int npx, int npy) {
    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;
    if ((i < npx) && (j < npy)) {
        int pos  = j * npx + i;
        int nvol = npx * npy;
        if (d_V[pos] > epsilon_h) {
            d_eta_max[pos] = max(d_eta_max[pos], d_V[pos] - d_V[3 * nvol + pos]);
        }
    }
}

__global__ void set_inicializacion(double* d_hmax, int npx, int npy) {
    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;
    if ((i < npx) && (j < npy)) {
        int pos = j * npx + i;
        d_hmax[pos] = -1e30;
    }
}

__global__ void set_inicializacion_bmin(double *d_bmin, int npx, int npy) {
    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;
    if ((i < npx) && (j < npy)) {
        int pos = j * npx + i;
        d_bmin[pos] = -1e30;
    }
}

__global__ void set_inicializacion_eta_max(double *d_eta_max, int npx, int npy) {
    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;
    if ((i < npx) && (j < npy)) {
        int pos = j * npx + i;
        d_eta_max[pos] = -1e30;
    }
}

static __device__ inline __attribute__((unused)) void spongeLayer(double *d_Vk, int nvol, int npx, int npy, double *d, int pos, int pos1, int pos2, int posd1, int posd2, double *CC, double *tCC, long unsigned int timeCC, double current_time) {

    double h1, h2, H1, H2;
    double2 q1, q2;

    // Corner
    if (pos2 > 0) {
        interpolate_timeSPL(CC, tCC, timeCC, current_time, h1, H1, q1, npx, npy, pos1);
        interpolate_timeSPL(CC, tCC, timeCC, current_time, h2, H2, q2, npx, npy, pos2);
        double aux = 0.5 * (d[posd1] + d[posd2]) * 0.5;
        d_Vk[pos]            = aux * (h1 + h2)     + (1.0 - d[posd1]) * d_Vk[pos] * 0.5            + (1.0 - d[posd2]) * d_Vk[pos] * 0.5;
        d_Vk[1 * nvol + pos] = aux * (q1.x + q2.x) + (1.0 - d[posd1]) * d_Vk[1 * nvol + pos] * 0.5 + (1.0 - d[posd2]) * d_Vk[1 * nvol + pos] * 0.5;
        d_Vk[2 * nvol + pos] = aux * (q1.y + q2.y) + (1.0 - d[posd1]) * d_Vk[2 * nvol + pos] * 0.5 + (1.0 - d[posd2]) * d_Vk[2 * nvol + pos] * 0.5;
        d_Vk[3 * nvol + pos] = aux * (H1 + H2)     + (1.0 - d[posd1]) * d_Vk[3 * nvol + pos] * 0.5 + (1.0 - d[posd2]) * d_Vk[3 * nvol + pos] * 0.5;
    }
    else {
        interpolate_timeSPL(CC, tCC, timeCC, current_time, h1, H1, q1, npx, npy, pos1);
        d_Vk[pos]            = d[posd1] * h1   + (1.0 - d[posd1]) * d_Vk[pos];
        d_Vk[1 * nvol + pos] = d[posd1] * q1.x + (1.0 - d[posd1]) * d_Vk[1 * nvol + pos];
        d_Vk[2 * nvol + pos] = d[posd1] * q1.y + (1.0 - d[posd1]) * d_Vk[2 * nvol + pos];
        d_Vk[3 * nvol + pos] = d[posd1] * H1   + (1.0 - d[posd1]) * d_Vk[3 * nvol + pos];
    }
}

// Untested!
static __device__ inline __attribute__((unused)) double dynamicRuptureInterpolation(double *dyrup, float *time_dyrup, double current_time, int pos, int npx, int npy) {

    float dt = time_dyrup[1] - time_dyrup[0];
    int t = int((current_time - time_dyrup[0]) / dt);
    t *= (t > 0);
    if (t > 40 - 1)
        return 0.0;
    int tp = t + 1;
    int posCC  = t * npx * npy + pos;
    int posCCp = tp * npx * npy + pos;

    double m = (dyrup[posCCp] - dyrup[posCC]) / (time_dyrup[tp] - time_dyrup[t]);
    double b = dyrup[posCC] - m * time_dyrup[t];
    return (m * current_time + b);
}

__global__ void states_updater(double* d_Vk, double* d_V0, double* d_V, int npx, int npy, int k, double epsilon_h, double invx, double invy, double *FL, double *FR, double *FU, double *FD, double *d, double current_time, double *d_CC, double *d_tCC, long unsigned int timeCC, double *dyrup, float *time_dyrup)
{
    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;

    if ((i < npx) && (j < npy)) {

        int pos  = j * npx + i;
        int nvol = npx * npy;

        // Planar flux updates across the 5 state planes
        d_Vk[pos]            -= invx * (FL[pos]            + FR[pos])            + invy * (FD[pos]            + FU[pos]);
        d_Vk[1 * nvol + pos] -= invx * (FL[1 * nvol + pos] + FR[1 * nvol + pos]) + invy * (FD[1 * nvol + pos] + FU[1 * nvol + pos]);
        d_Vk[2 * nvol + pos] -= invx * (FL[2 * nvol + pos] + FR[2 * nvol + pos]) + invy * (FD[2 * nvol + pos] + FU[2 * nvol + pos]);
        d_Vk[3 * nvol + pos] -= invx * (FL[3 * nvol + pos] + FR[3 * nvol + pos]) + invy * (FD[3 * nvol + pos] + FU[3 * nvol + pos]);
        d_Vk[4 * nvol + pos] -= invx * (FL[4 * nvol + pos] + FR[4 * nvol + pos]) + invy * (FD[4 * nvol + pos] + FU[4 * nvol + pos]);

#if SPL > 0
        if (current_time < d_tCC[timeCC - 2]) {
            // Strip down
            if (j < SPL) {
                // bottom left corner
                if (i < SPL) {
                    spongeLayer(d_Vk, nvol, npx, npy, d, pos, i, 2 * npx + npy + j, i, j, d_CC, d_tCC, timeCC, current_time);
                }
                // Bottom right corner
                else if (i > npx - SPL) {
                    spongeLayer(d_Vk, nvol, npx, npy, d, pos, i, npx + j, npx - 1 - i, j, d_CC, d_tCC, timeCC, current_time);
                }
                // Central strip
                else {
                    spongeLayer(d_Vk, nvol, npx, npy, d, pos, i, -1, j, -1, d_CC, d_tCC, timeCC, current_time);
                }
            }
            else {
                // Upper strip
                if (j > npy - SPL) {
                    // Up left corner
                    if (i < SPL) {
                        spongeLayer(d_Vk, nvol, npx, npy, d, pos, npx + npy + i, 2 * npx + npy + j, i, npy - 1 - j, d_CC, d_tCC, timeCC, current_time);
                    }
                    // Up right corner
                    else if (i > npx - SPL) {
                        spongeLayer(d_Vk, nvol, npx, npy, d, pos, npx + npy + i, npx + j, npx - 1 - i, npy - 1 - j, d_CC, d_tCC, timeCC, current_time);
                    }
                    // Central strip
                    else {
                        spongeLayer(d_Vk, nvol, npx, npy, d, pos, npx + npy + i, -1, npy - 1 - j, -1, d_CC, d_tCC, timeCC, current_time);
                    }
                }
                else {
                    // Left strip
                    if (i < SPL) {
                        spongeLayer(d_Vk, nvol, npx, npy, d, pos, 2 * npx + npy + j, -1, i, -1, d_CC, d_tCC, timeCC, current_time);
                    }
                    // Right strip
                    if (i > npx - SPL) {
                        spongeLayer(d_Vk, nvol, npx, npy, d, pos, npx + j, -1, npx - 1 - i, -1, d_CC, d_tCC, timeCC, current_time);
                    }
                }
            }
        }
#endif

#if RK == 1
        double coef0 = 0.0 * (k == 0) + 0.5 * (k == 1);
        double coef1 = 1.0 - coef0;

        double2 q;
        double h   = coef0 * d_V0[pos]            + coef1 * d_Vk[pos];
        q.x        = coef0 * d_V0[1 * nvol + pos] + coef1 * d_Vk[1 * nvol + pos];
        q.y        = coef0 * d_V0[2 * nvol + pos] + coef1 * d_Vk[2 * nvol + pos];
        double z   = coef0 * d_V0[3 * nvol + pos] + coef1 * d_Vk[3 * nvol + pos];
        double hc  = coef0 * d_V0[4 * nvol + pos] + coef1 * d_Vk[4 * nvol + pos];

        state_filter(h, q, hc, epsilon_h);

        d_Vk[pos]            = h;
        d_Vk[1 * nvol + pos] = q.x;
        d_Vk[2 * nvol + pos] = q.y;
        d_Vk[3 * nvol + pos] = z;
        d_Vk[4 * nvol + pos] = hc;

        if (k == order - 1) {
            d_V0[pos]            = h;
            d_V0[1 * nvol + pos] = q.x;
            d_V0[2 * nvol + pos] = q.y;
            d_V0[3 * nvol + pos] = z;
            d_V0[4 * nvol + pos] = hc;

            d_V[pos]             = h;
            d_V[1 * nvol + pos]  = q.x;
            d_V[2 * nvol + pos]  = q.y;
            d_V[3 * nvol + pos]  = z;
            d_V[4 * nvol + pos]  = hc;
        }
        else {
            d_V[pos]             = h;
            d_V[1 * nvol + pos]  = q.x;
            d_V[2 * nvol + pos]  = q.y;
            d_V[3 * nvol + pos]  = z;
            d_V[4 * nvol + pos]  = hc;
        }
#else
        double h,z,hc;
        double2 q;
        h = d_Vk[pos];
        q.x = d_Vk[1 * nvol + pos];
        q.y = d_Vk[2 * nvol + pos];
        z = d_Vk[3 * nvol + pos];
        hc = d_Vk[4 * nvol + pos];

        state_filter(h, q, hc, epsilon_h);

        d_V[pos]            = d_V0[pos]            = h;
        d_V[1 * nvol + pos] = d_V0[1 * nvol + pos] = q.x;
        d_V[2 * nvol + pos] = d_V0[2 * nvol + pos] = q.y;
        d_V[3 * nvol + pos] = d_V0[3 * nvol + pos] = d_Vk[3 * nvol + pos];
        d_V[4 * nvol + pos] = d_V0[4 * nvol + pos] = hc;
#endif
    }
}
