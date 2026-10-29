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

static __device__ inline __attribute__((unused)) double avg(double a, double b)
{
    double d, aa, bb;
    aa = fabs(a);
    bb = fabs(b);
    if ((aa + bb) > EPSILON) {
        d = (a * bb + b * aa) / (aa + bb);
    }
    else {
        d = 0.0;
    }
    return d;
}

static __device__ inline void friction(double* d_Vk, double* d_V0, int pos, int nvol, double epsilon_h, double bottom_friction, double dt) {

    double h;
    double2 q, u;

    h   = d_V0[pos];
    q.x = d_V0[1 * nvol + pos];
    q.y = d_V0[2 * nvol + pos];

    u_calculator(h, q, epsilon_h, u);
    double modu = sqrt(u.x * u.x + u.y * u.y);

    double c = dt * bottom_friction * modu / (h + EPSILON);

#if Manning == 1
    c = dt * bottom_friction * bottom_friction * modu / (pow(h, 4.0 / 3.0));
    d_Vk[1 * nvol + pos] *= c;
    d_Vk[2 * nvol + pos] *= c;
#else
    d_Vk[1 * nvol + pos] /= (1.0 + c);
    d_Vk[2 * nvol + pos] /= (1.0 + c);
#endif
}

static __device__ inline void erosion_deposition(double* d_Vk, double* d_V0, const SimParams &params, double dt, int pos){

    double Fe, Fd, phi;
    double h, c;
    double vs, Rp, nu;
    double2 q, u;

    nu = 1.002e-6;          // Water kinematic viscosity
    double s = params.R - 1.0; // Difference of relative densities
    int nvol = params.nvol;

    h   = d_V0[pos];
    q.x = d_V0[1 * nvol + pos];
    q.y = d_V0[2 * nvol + pos];
    c   = d_V0[4 * nvol + pos] / (h + params.epsilon_h);

    u_calculator(h, q, params.epsilon_h, u);
    double modu = sqrt(u.x * u.x + u.y * u.y);

    vs = sqrt((13.95 * nu / params.DS) * (13.95 * nu / params.DS) + 1.09 * s * G * params.DS) - 13.95 * nu / params.DS; // Vertical vel
    Rp = (params.DS * sqrt(s * G * params.DS)) / nu; // Particle Reynolds number

    calculate_erosion(u, modu, Rp, vs, params, Fe);
    calculate_deposition(c, vs, Rp, params, Fd);

    phi = (Fe - Fd) * SUSPENSION;



    d_Vk[3 * nvol + pos] -= dt * phi;
    d_Vk[4 * nvol + pos] += dt * phi;
}

__global__ void reconstructor(double* d_V, double* reconsL, double* reconsR, double* reconsD, double* reconsU, double delta_T, double inv_dx, double inv_dy, SimParams params, double* d_Vk, double *d_V0) {

    int i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    int j = blockIdx.y * NUM_THREADS_VERTICAL_VOL   + threadIdx.y;

    if ((i < params.npx) && (j < params.npy)) {
        int pos  = j * params.npx + i;
        int nvol = params.nvol;

#if order == 1


        reconsL[pos]            = reconsR[pos]            = reconsU[pos]            = reconsD[pos]            = d_V[pos];
        reconsL[1 * nvol + pos] = reconsR[1 * nvol + pos] = reconsU[1 * nvol + pos] = reconsD[1 * nvol + pos] = d_V[1 * nvol + pos];
        reconsL[2 * nvol + pos] = reconsR[2 * nvol + pos] = reconsU[2 * nvol + pos] = reconsD[2 * nvol + pos] = d_V[2 * nvol + pos];
        reconsL[3 * nvol + pos] = reconsR[3 * nvol + pos] = reconsU[3 * nvol + pos] = reconsD[3 * nvol + pos] = d_V[3 * nvol + pos];
        reconsL[4 * nvol + pos] = reconsR[4 * nvol + pos] = reconsU[4 * nvol + pos] = reconsD[4 * nvol + pos] = d_V[4 * nvol + pos];


#else

        int posleft, posright;
        int im, ip, jm, jp;
        double dta, Hl, Hi, Hr;
        double hi, hl, hr;
        double2 ui, ul, ur;
        double dl, dr, du;
        double ci, cr, cl;
        double deta;

        // **********************
        // x direction
        // **********************

        hi = d_V[pos];
        Hi = d_V[3 * nvol + pos];
        ci = d_V[4 * nvol + pos] / (hi + EPSILON);

        im = i - 1 * (i > 0);
        ip = i + 1 * (i < params.npx - 1);
#if WBC == 2
        if (i == 0) im = params.npx - 1;
        if (i == params.npx - 1) ip = 0;
#endif
        posleft  = j * params.npx + im;
        posright = j * params.npx + ip;

        hl = d_V[posleft];
        Hl = d_V[3 * nvol + posleft];
        cl = d_V[4 * nvol + posleft] / (hl + EPSILON);

        hr = d_V[posright];
        Hr = d_V[3 * nvol + posright];
        cr = d_V[4 * nvol + posright] / (hr + EPSILON);
        dta = delta_T * inv_dx;

        // h reconstruction
        dl = hi - hl;
        dr = hr - hi;
        double dh = avg(dl, dr);
        if (0.5 * fabs(dh) > hi) {
            dh = 2.0 * hi * SGN(dh);
        }
        reconsL[pos] = hi - 0.5 * dh;
        reconsR[pos] = hi + 0.5 * dh;

        // eta and zb reconstruction
        double zm;
        zm = min(Hl, Hi);
        dl = max(hi - Hi + zm, 0.0) - max(hl - Hl + zm, 0.0);
        zm = min(Hr, Hi);
        dr = max(hr - Hr + zm, 0.0) - max(hi - Hi + zm, 0.0);
        deta = avg(dl, dr);
        double dz = dh - deta;

        reconsL[3 * nvol + pos] = Hi - 0.5 * dz;
        reconsR[3 * nvol + pos] = Hi + 0.5 * dz;

        // u reconstruction
        u_calculator_double(hi, d_V[1 * nvol + pos], d_V[2 * nvol + pos], params.epsilon_h, ui);
        u_calculator_double(hl, d_V[1 * nvol + posleft], d_V[2 * nvol + posleft], params.epsilon_h, ul);
        u_calculator_double(hr, d_V[1 * nvol + posright], d_V[2 * nvol + posright], params.epsilon_h, ur);
        dl = ui.x - ul.x;
        dr = ur.x - ui.x;
        du = avg(dl, dr);
        double dq = hi * du + ui.x * dh;

        reconsL[1 * nvol + pos] = d_V[1 * nvol + pos] - 0.5 * dq;
        reconsR[1 * nvol + pos] = d_V[1 * nvol + pos] + 0.5 * dq;

        dl = ui.y - ul.y;
        dr = ur.y - ui.y;
        du = avg(dl, dr);
        dq = hi * du + ui.y * dh;

        reconsL[2 * nvol + pos] = d_V[2 * nvol + pos] - 0.5 * dq;
        reconsR[2 * nvol + pos] = d_V[2 * nvol + pos] + 0.5 * dq;

        // hc reconstruction
        dl = ci - cl;
        dr = cr - ci;
        du = avg(dl, dr);

        reconsL[4 * nvol + pos] = (ci - 0.5 * du) * reconsL[pos];
        reconsR[4 * nvol + pos] = (ci + 0.5 * du) * reconsR[pos];

        // Integral
        d_Vk[1 * nvol + pos] -= dta * G * hi * deta;

#if RK == 0
        double Tpx = G * hi * deta;
#endif

        // **********************
        // y direction
        // **********************

        jm = j - 1 * (j > 0);
        jp = j + 1 * (j < params.npy - 1);
#if WBC == 2
        if (j == 0) jm = params.npy - 1;
        if (j == params.npy - 1) jp = 0;
#endif
        posleft  = jm * params.npx + i;
        posright = jp * params.npx + i;

        hl = d_V[posleft];
        Hl = d_V[3 * nvol + posleft];
        cl = d_V[4 * nvol + posleft] / (hl + EPSILON);

        hr = d_V[posright];
        Hr = d_V[3 * nvol + posright];
        cr = d_V[4 * nvol + posright] / (hr + EPSILON);
        dta = delta_T * inv_dy;

        // h reconstruction
        dl = hi - hl;
        dr = hr - hi;
        dh = avg(dl, dr);
        if (0.5 * fabs(dh) > hi) {
            dh = 2.0 * hi * SGN(dh);
        }
        reconsD[pos] = hi - 0.5 * dh;
        reconsU[pos] = hi + 0.5 * dh;

        // eta and zb reconstruction
        zm = min(Hl, Hi);
        dl = max(hi - Hi + zm, 0.0) - max(hl - Hl + zm, 0.0);
        zm = min(Hr, Hi);
        dr = max(hr - Hr + zm, 0.0) - max(hi - Hi + zm, 0.0);
        deta = avg(dl, dr);
        dz = dh - deta;

        reconsD[3 * nvol + pos] = Hi - 0.5 * dz;
        reconsU[3 * nvol + pos] = Hi + 0.5 * dz;

        // u reconstruction
        u_calculator_double(hl, d_V[1 * nvol + posleft], d_V[2 * nvol + posleft], params.epsilon_h, ul);
        u_calculator_double(hr, d_V[1 * nvol + posright], d_V[2 * nvol + posright], params.epsilon_h, ur);
        dl = ui.x - ul.x;
        dr = ur.x - ui.x;
        du = avg(dl, dr);
        dq = hi * du + ui.x * dh;

        reconsD[1 * nvol + pos] = d_V[1 * nvol + pos] - 0.5 * dq;
        reconsU[1 * nvol + pos] = d_V[1 * nvol + pos] + 0.5 * dq;

        dl = ui.y - ul.y;
        dr = ur.y - ui.y;
        du = avg(dl, dr);
        dq = hi * du + ui.y * dh;

        reconsD[2 * nvol + pos] = d_V[2 * nvol + pos] - 0.5 * dq;
        reconsU[2 * nvol + pos] = d_V[2 * nvol + pos] + 0.5 * dq;

        // hc reconstruction
        dl = ci - cl;
        dr = cr - ci;
        du = avg(dl, dr);

        reconsD[4 * nvol + pos] = (ci - 0.5 * du) * reconsD[pos];
        reconsU[4 * nvol + pos] = (ci + 0.5 * du) * reconsU[pos];

        // Integral
        d_Vk[2 * nvol + pos] -= dta * G * hi * deta;

#if RK == 0
        double Tpy = G * hi * deta;
#endif

#if RK == 0
        double2 uRR, uL, uU, uD;
        u_calculator_double(reconsR[pos], reconsR[1 * nvol + pos], reconsR[2 * nvol + pos], params.epsilon_h, uRR);
        u_calculator_double(reconsL[pos], reconsL[1 * nvol + pos], reconsL[2 * nvol + pos], params.epsilon_h, uL);
        u_calculator_double(reconsU[pos], reconsU[1 * nvol + pos], reconsU[2 * nvol + pos], params.epsilon_h, uU);
        u_calculator_double(reconsD[pos], reconsD[1 * nvol + pos], reconsD[2 * nvol + pos], params.epsilon_h, uD);

        double qhx, qhy;
        qhx = (reconsR[1 * nvol + pos] - reconsL[1 * nvol + pos]) * inv_dx;
        qhy = (reconsU[2 * nvol + pos] - reconsD[2 * nvol + pos]) * inv_dy;

        // MUSCL-Hancock half-step predictor
        double update_h = 0.5 * delta_T * (qhx + qhy);
        reconsL[pos] -= update_h;
        reconsR[pos] -= update_h;
        reconsD[pos] -= update_h;
        reconsU[pos] -= update_h;

        qhx = (reconsR[1 * nvol + pos] * uRR.x - reconsL[1 * nvol + pos] * uL.x) * inv_dx;
        qhy = (reconsU[1 * nvol + pos] * uU.y  - reconsD[1 * nvol + pos] * uD.y)  * inv_dy;
        double update_qx = 0.5 * delta_T * (qhx + qhy);
        reconsL[1 * nvol + pos] -= update_qx;
        reconsR[1 * nvol + pos] -= update_qx;
        reconsD[1 * nvol + pos] -= update_qx;
        reconsU[1 * nvol + pos] -= update_qx;

        qhx = (reconsR[1 * nvol + pos] * uRR.y - reconsL[1 * nvol + pos] * uL.y) * inv_dx;
        qhy = (reconsU[2 * nvol + pos] * uU.y  - reconsD[2 * nvol + pos] * uD.y)  * inv_dy;
        double update_qy = 0.5 * delta_T * (qhx + qhy);
        reconsL[2 * nvol + pos] -= update_qy;
        reconsR[2 * nvol + pos] -= update_qy;
        reconsD[2 * nvol + pos] -= update_qy;
        reconsU[2 * nvol + pos] -= update_qy;

        qhx = (bedload(reconsR[pos], uRR, 0.0, params) - bedload(reconsL[pos], uL, 0.0, params)) * inv_dx;
        qhy = (bedload(reconsU[pos], uU, 0.0, params)  - bedload(reconsD[pos], uD, 0.0, params)) * inv_dy;
        double update_zb = 0.5 * delta_T * (qhx + qhy);
        reconsL[3 * nvol + pos] -= update_zb;
        reconsR[3 * nvol + pos] -= update_zb;
        reconsD[3 * nvol + pos] -= update_zb;
        reconsU[3 * nvol + pos] -= update_zb;

        qhx = (reconsR[4 * nvol + pos] * uRR.x - reconsL[4 * nvol + pos] * uL.x) * inv_dx;
        qhy = (reconsU[4 * nvol + pos] * uU.y  - reconsD[4 * nvol + pos] * uD.y)  * inv_dy;
        double update_hc = 0.5 * delta_T * (qhx + qhy);
        reconsL[4 * nvol + pos] -= update_hc;
        reconsR[4 * nvol + pos] -= update_hc;
        reconsD[4 * nvol + pos] -= update_hc;
        reconsU[4 * nvol + pos] -= update_hc;

        double update_Tpx = 0.5 * delta_T * inv_dx * Tpx;
        reconsL[1 * nvol + pos] -= update_Tpx;
        reconsU[1 * nvol + pos] -= update_Tpx;
        reconsD[1 * nvol + pos] -= update_Tpx;
        reconsR[1 * nvol + pos] -= update_Tpx;

        double update_Tpy = 0.5 * delta_T * inv_dy * Tpy;
        reconsL[2 * nvol + pos] -= update_Tpy;
        reconsU[2 * nvol + pos] -= update_Tpy;
        reconsD[2 * nvol + pos] -= update_Tpy;
        reconsR[2 * nvol + pos] -= update_Tpy;

#endif

        // Filters
        state_filter_double(reconsL[pos], reconsL[1 * nvol + pos], reconsL[2 * nvol + pos], reconsL[4 * nvol + pos], params.epsilon_h);
        state_filter_double(reconsR[pos], reconsR[1 * nvol + pos], reconsR[2 * nvol + pos], reconsR[4 * nvol + pos], params.epsilon_h);
        state_filter_double(reconsU[pos], reconsU[1 * nvol + pos], reconsU[2 * nvol + pos], reconsU[4 * nvol + pos], params.epsilon_h);
        state_filter_double(reconsD[pos], reconsD[1 * nvol + pos], reconsD[2 * nvol + pos], reconsD[4 * nvol + pos], params.epsilon_h);

#endif

        // Frictions and source terms
        friction(d_Vk, d_V0, pos, nvol, params.epsilon_h, params.FR, delta_T);
        erosion_deposition(d_Vk, d_V0, params, delta_T, pos);
    }
}
