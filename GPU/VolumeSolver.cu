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

static __device__ __forceinline__ void RiemannSolver_Vols(double hhl, double zb0, double2 qql, double hcl, double hhr, double zb1, double2 qqr, double hcr, int pos0, int pos1, const SimParams &params, double *FL, double *FR, int normal_x, int normal_y, int fixed_sediment) {

    double2 q0, q1;
    double2 ql, qr;
    double2 u0, u1;
    double h0, h1;
    double c0, c1, hc0, hc1;

    // x = normal, y = tangential
    ql.x = qql.x * normal_x + qql.y * normal_y;
    ql.y = qql.y * normal_x - qql.x * normal_y;
    qr.x = qqr.x * normal_x + qqr.y * normal_y;
    qr.y = qqr.y * normal_x - qqr.x * normal_y;
    u_calculator(hhl, ql, params.epsilon_h, u0);
    u_calculator(hhr, qr, params.epsilon_h, u1);
    c0 = hcl / (hhl + EPSILON);
    c1 = hcr / (hhr + EPSILON);

    // Hydrostatic Reconstruction
    double zmH = min(zb0, zb1);
    h0 = max(hhl - zb0 + zmH, 0.0);
    h1 = max(hhr - zb1 + zmH, 0.0);

    q0.x = u0.x * h0;
    q0.y = u0.y * h0;

    q1.x = u1.x * h1;
    q1.y = u1.y * h1;

    hc0 = h0 * c0;
    hc1 = h1 * c1;

    double deta = h1 - h0;

    double F0[5], F1[5];
    F0[0] = F0[1] = F0[2] = F0[3] = F0[4] = 0.0;
    F1[0] = F1[1] = F1[2] = F1[3] = F1[4] = 0.0;
    Cflux(h0, u0, hc0, F0, fixed_sediment, zb1 - zb0, params);
    Cflux(h1, u1, hc1, F1, fixed_sediment, zb1 - zb0, params);

    double hm, zm; // Roe average
    double2 um;
    roe_average(h0, u0, zb0, h1, u1, zb1, hm, um, zm);

    double eigenvalues[3], SL, SR;

    calc_eigen_aprox0(h0, u0, eigenvalues);
    SL = eigenvalues[0];

    calc_eigen_aprox0(h1, u1, eigenvalues);
    SR = eigenvalues[2];

    calc_eigen_aprox0(hm, um, eigenvalues);
    SL = min(SL, eigenvalues[0]);
    SR = max(SR, eigenvalues[2]);

    double a0, a1;
    if (SR - SL > EPSILON) {
        a0 = (SR * fabs(SL) - SL * fabs(SR)) / (SR - SL);
        a1 = (fabs(SR) - fabs(SL)) / (SR - SL);
    }
    else {
        a0 = max(fabs(SL), fabs(SR));
        a1 = 0.0;
    }

    double df[5], vis[5], DU[5];
    double s0, s1, s2;

    for (int i = 0; i < 5; i++) {
        df[i] = F1[i] - F0[i];
    }

    df[1] += G * hm * deta;

    DU[0] = h1 - h0;
    DU[1] = q1.x - q0.x;
    DU[2] = q1.y - q0.y;
    DU[3] = zb1 - zb0;
    DU[4] = hc1 - hc0;

    s0 = SGN(zb1 - zb0);
    s1 = SGN(df[3]);
    s2 = SGN(a1);

    if (fabs(s0) > EPSILON) {
        s1 = s1 / s0;
    }
    else {
        s1 = 0.0;
    }

    for (int i = 0; i < 5; i++) {
        vis[i] = a0 * DU[i] + a1 * df[i];
    }

    if (s2 == s1) {
        vis[3] = a1 * df[3];
    }
    else {
        vis[3] = -a1 * df[3];
    }

    double Dm[5], Dp[5];
    for (int i = 0; i < 5; i++) {
        Dm[i] = 0.5 * (df[i] - vis[i]) + F0[i];
        Dp[i] = 0.5 * (df[i] + vis[i]) - F1[i];
    }
    int nvol = params.nvol;
#if WBC == 2
    FR[pos0]            = Dm[0];
    FR[1 * nvol + pos0] = (normal_x * Dm[1] - normal_y * Dm[2]);
    FR[2 * nvol + pos0] = (normal_y * Dm[1] + normal_x * Dm[2]);
    FR[3 * nvol + pos0] = Dm[3];
    FR[4 * nvol + pos0] = Dm[4];

    FL[pos1]            = Dp[0];
    FL[1 * nvol + pos1] = (normal_x * Dp[1] - normal_y * Dp[2]);
    FL[2 * nvol + pos1] = (normal_y * Dp[1] + normal_x * Dp[2]);
    FL[3 * nvol + pos1] = Dp[3];
    FL[4 * nvol + pos1] = Dp[4];
#else
    if (pos1 >= 0) {
        FR[pos0]            = Dm[0];
        FR[1 * nvol + pos0] = (normal_x * Dm[1] - normal_y * Dm[2]);
        FR[2 * nvol + pos0] = (normal_y * Dm[1] + normal_x * Dm[2]);
        FR[3 * nvol + pos0] = Dm[3];
        FR[4 * nvol + pos0] = Dm[4];

        FL[pos1]            = Dp[0];
        FL[1 * nvol + pos1] = (normal_x * Dp[1] - normal_y * Dp[2]);
        FL[2 * nvol + pos1] = (normal_y * Dp[1] + normal_x * Dp[2]);
        FL[3 * nvol + pos1] = Dp[3];
        FL[4 * nvol + pos1] = Dp[4];
    }
    else {
        if (pos1 == -1) {
            FL[pos0]            = Dm[0];
            FL[1 * nvol + pos0] = (normal_x * Dm[1] - normal_y * Dm[2]);
            FL[2 * nvol + pos0] = (normal_y * Dm[1] + normal_x * Dm[2]);
            FL[3 * nvol + pos0] = Dm[3];
            FL[4 * nvol + pos0] = Dm[4];
        }
        else {
            FR[pos0]            = Dm[0];
            FR[1 * nvol + pos0] = (normal_x * Dm[1] - normal_y * Dm[2]);
            FR[2 * nvol + pos0] = (normal_y * Dm[1] + normal_x * Dm[2]);
            FR[3 * nvol + pos0] = Dm[3];
            FR[4 * nvol + pos0] = Dm[4];
        }
    }
#endif
}
static __device__ inline __attribute__((unused)) void interpolate_time(double *CC, double *tCC, long unsigned int timeCC, double current_time, double &h1, double &zb1, double2 &q1, double &hc1, int npx, int npy, int pos, double *tidesBd, float *tidesT, int tidesNPT) {

    double dt = tCC[1] - tCC[0];

    int t = int((current_time - tCC[0]) / dt);
    t *= (t > 0);
    int tp = t + 1;
    int posCC  = t * (2 * npx + 2 * npy) * 4 + pos * 4;
    int posCCp = tp * (2 * npx + 2 * npy) * 4 + pos * 4;

    double m = (CC[posCCp] - CC[posCC]) / dt;
    double b = CC[posCC] - m * tCC[t];
    h1 = m * current_time + b;

    m = (CC[posCCp + 3] - CC[posCC + 3]) / dt;
    b = CC[posCC + 3] - m * tCC[t];
    zb1 = m * current_time + b;

    h1 += zb1;

    m = (CC[posCCp + 1] - CC[posCC + 1]) / dt;
    b = CC[posCC + 1] - m * tCC[t];
    q1.x = m * current_time + b;
    q1.x *= h1;

    m = (CC[posCCp + 2] - CC[posCC + 2]) / dt;
    b = CC[posCC + 2] - m * tCC[t];
    q1.y = m * current_time + b;
    q1.y *= h1;

    hc1 = 0.0;

    // tides
    // dt = tidesT[1]-tidesT[0];
    // t = int((current_time-tidesT[0])/dt) ;
    // t*=(t>0);
    // tp = t+1;
    // posCC = t*(2*npx+2*npy) + pos;
    // posCCp = tp*(2*npx+2*npy) + pos;
    //
    // m = (tidesBd[posCCp]-tidesBd[posCC])/dt;
    // b = tidesBd[posCC] - m*tidesT[t];
    // h1 = h1 + (m*current_time+b) * (current_time<tidesT[tidesNPT-1]);
}

__global__ void VolumeSolver(int *fixH, double current_time, SimParams params, double* reconsL, double* reconsR, double* reconsD, double* reconsU, double *FL, double *FR, double *FU, double *FD, double *d_CC, double *d_tCC, long unsigned int timeCC, double *tidesBd, float *tidesT, int tidesNPT) {

    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;
    int pos0, pos1;
    int nvol = params.nvol;
    double h0, h1, zb0, zb1, hc0, hc1;
    double2 q0, q1;
    int fixed_sediment = 0;

    if ((i <= params.npx) && (j <= params.npy)) {
        if ((i <= params.npx) && (j < params.npy)) {
            // Riemann problems at vertical edges
            pos0 = i - 1 + j * params.npx;
            pos1 = i + j * params.npx;

            if (i == 0) {
#if WBC == 2

#else
                pos0 = pos1;
                pos1 = -1;
                h1   = h0   = reconsL[pos0];
                q1.x = q0.x = reconsL[1 * nvol + pos0];
                q1.y = q0.y = reconsL[2 * nvol + pos0];
                zb1  = zb0  = reconsL[3 * nvol + pos0];
                hc1  = hc0  = reconsL[4 * nvol + pos0];

#if WBC == 1
                q1.x = -q0.x;
#endif
#if (SPL <= 0 && Read_files == 1)
                if (current_time < d_tCC[timeCC - 2]) {
                    interpolate_time(d_CC, d_tCC, timeCC, current_time, h1, zb1, q1, hc1, params.npx, params.npy, 2 * params.npx + params.npy + j, tidesBd, tidesT, tidesNPT);
                }
                else {
                    q1.x = -q0.x;
                }
#endif
                fixed_sediment = fixH[pos0];
                RiemannSolver_Vols(h0, zb0, q0, hc0, h1, zb1, q1, hc1, pos0, pos1, params, FL, FR, -1, 0, fixed_sediment);
#endif
            }
            else {
                if (i == params.npx) {
#if WBC == 2
                    pos1 = j * params.npx;

                    h0   = reconsR[pos0];
                    q0.x = reconsR[1 * nvol + pos0];
                    q0.y = reconsR[2 * nvol + pos0];
                    zb0  = reconsR[3 * nvol + pos0];
                    hc0  = reconsR[4 * nvol + pos0];

                    h1   = reconsL[pos1];
                    q1.x = reconsL[1 * nvol + pos1];
                    q1.y = reconsL[2 * nvol + pos1];
                    zb1  = reconsL[3 * nvol + pos1];
                    hc1  = reconsL[4 * nvol + pos1];
#else
                    pos1 = -2;

                    h1   = h0   = reconsR[pos0];
                    q1.x = q0.x = reconsR[1 * nvol + pos0];
                    q1.y = q0.y = reconsR[2 * nvol + pos0];
                    zb1  = zb0  = reconsR[3 * nvol + pos0];
                    hc1  = hc0  = reconsR[4 * nvol + pos0];

#if WBC == 1
                    q1.x = -q0.x;
#endif
#endif
#if (SPL <= 0 && Read_files == 1)
                    if (current_time < d_tCC[timeCC - 2]) {
                        interpolate_time(d_CC, d_tCC, timeCC, current_time, h1, zb1, q1, hc1, params.npx, params.npy, params.npx + j, tidesBd, tidesT, tidesNPT);
                    }
                    else {
                        q1.x = -q0.x;
                    }
#endif
                    fixed_sediment = fixH[pos0];
                    RiemannSolver_Vols(h0, zb0, q0, hc0, h1, zb1, q1, hc1, pos0, pos1, params, FL, FR, 1, 0, fixed_sediment);
                }
                else {
                    h0   = reconsR[pos0];
                    q0.x = reconsR[1 * nvol + pos0];
                    q0.y = reconsR[2 * nvol + pos0];
                    zb0  = reconsR[3 * nvol + pos0];
                    hc0  = reconsR[4 * nvol + pos0];

                    h1   = reconsL[pos1];
                    q1.x = reconsL[1 * nvol + pos1];
                    q1.y = reconsL[2 * nvol + pos1];
                    zb1  = reconsL[3 * nvol + pos1];
                    hc1  = reconsL[4 * nvol + pos1];

                    fixed_sediment = fixH[pos0] * fixH[pos1];
                    RiemannSolver_Vols(h0, zb0, q0, hc0, h1, zb1, q1, hc1, pos0, pos1, params, FL, FR, 1, 0, fixed_sediment);
                }
            }
        }

        if ((i < params.npx) && (j <= params.npy)) {
            // Riemann problems at horizontal edges
            pos0 = i + (j - 1) * params.npx;
            pos1 = i + j * params.npx;

            if (j == 0) {
#if WBC == 2

#else
                pos0 = pos1;
                pos1 = -1;

                h1   = h0   = reconsD[pos0];
                q1.x = q0.x = reconsD[1 * nvol + pos0];
                q1.y = q0.y = reconsD[2 * nvol + pos0];
                zb1  = zb0  = reconsD[3 * nvol + pos0];
                hc1  = hc0  = reconsD[4 * nvol + pos0];
#if WBC == 1
                q1.y = -q0.y;
#endif
#if (SPL <= 0 && Read_files == 1)
                if (current_time < d_tCC[timeCC - 2]) {
                    interpolate_time(d_CC, d_tCC, timeCC, current_time, h1, zb1, q1, hc1, params.npx, params.npy, i, tidesBd, tidesT, tidesNPT);
                }
                else {
                    q1.y = -q0.y;
                }
#endif
                fixed_sediment = fixH[pos0];
                RiemannSolver_Vols(h0, zb0, q0, hc0, h1, zb1, q1, hc1, pos0, pos1, params, FD, FU, 0, -1, fixed_sediment);
#endif
            }
            else {
                if (j == params.npy) {
#if WBC == 2
                    pos1 = i;

                    h0   = reconsU[pos0];
                    q0.x = reconsU[1 * nvol + pos0];
                    q0.y = reconsU[2 * nvol + pos0];
                    zb0  = reconsU[3 * nvol + pos0];
                    hc0  = reconsU[4 * nvol + pos0];

                    h1   = reconsD[pos1];
                    q1.x = reconsD[1 * nvol + pos1];
                    q1.y = reconsD[2 * nvol + pos1];
                    zb1  = reconsD[3 * nvol + pos1];
                    hc1  = reconsD[4 * nvol + pos1];
#else
                    pos1 = -2;
                    h1   = h0   = reconsU[pos0];
                    q1.x = q0.x = reconsU[1 * nvol + pos0];
                    q1.y = q0.y = reconsU[2 * nvol + pos0];
                    zb1  = zb0  = reconsU[3 * nvol + pos0];
                    hc1  = hc0  = reconsU[4 * nvol + pos0];
#if WBC == 1
                    q1.y = -q0.y;
#endif
#endif
#if (SPL <= 0 && Read_files == 1)
                    if (current_time < d_tCC[timeCC - 2]) {
                        interpolate_time(d_CC, d_tCC, timeCC, current_time, h1, zb1, q1, hc1, params.npx, params.npy, params.npx + params.npy + i, tidesBd, tidesT, tidesNPT);
                    }
                    else {
                        q1.y = -q0.y;
                    }
#endif
                    fixed_sediment = fixH[pos0];
                    RiemannSolver_Vols(h0, zb0, q0, hc0, h1, zb1, q1, hc1, pos0, pos1, params, FD, FU, 0, 1, fixed_sediment);
                }
                else {
                    h0   = reconsU[pos0];
                    q0.x = reconsU[1 * nvol + pos0];
                    q0.y = reconsU[2 * nvol + pos0];
                    zb0  = reconsU[3 * nvol + pos0];
                    hc0  = reconsU[4 * nvol + pos0];

                    h1   = reconsD[pos1];
                    q1.x = reconsD[1 * nvol + pos1];
                    q1.y = reconsD[2 * nvol + pos1];
                    zb1  = reconsD[3 * nvol + pos1];
                    hc1  = reconsD[4 * nvol + pos1];

                    fixed_sediment = fixH[pos0] * fixH[pos1];
                    RiemannSolver_Vols(h0, zb0, q0, hc0, h1, zb1, q1, hc1, pos0, pos1, params, FD, FU, 0, 1, fixed_sediment);
                }
            }
        }
    }
}
