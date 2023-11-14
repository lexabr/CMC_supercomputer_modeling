#include <iostream>
#include <cmath>
#include <vector>

int N, Nx, Ny, Nz;
double Lx, Ly, Lz;
double hx, hy, hz;

const double a_2 = 1/9;
const double tau = 0.001;
const int K = 20;
const int num_stages = 3;


double func_u(double x, double y, double z, double t) {
    double a_t = M_PI / 3 * sqrt(4 / pow(Lx, 2) + 1 / pow(Ly, 2) + 4 / pow(Lz, 2));

    return sin(2 * M_PI * x / Lx) * sin(M_PI * y / Ly + M_PI) * sin(2 * M_PI * z / Lz + 2 * M_PI) * cos(a_t * t + M_PI);
}

double func_phi(double x, double y, double z) {
    return func_u(x, y, z, 0);
}

double real_x(int i) {return i * hx;}
double real_y(int j) {return j * hy;}
double real_z(int k) {return k * hz;}
double real_t(int t) {return t * tau;}


int index(int i, int j, int k) {
    return Ny * Nx * k + Nx * j + i;
}

auto idx_to_ijk(int idx) {
    int i = idx % Nx;
    int j = (idx / Nx) % Ny;
    int k = (idx / Nx / Ny) % Nz;

    return std::make_tuple(i, j, k);
}


double laplace(const double **data, int st, int i, int j, int k) {
    
    return data[st][index(i - 1, j, k)] - 2 * data[st][index(i, j, k)] + data[st][index(i + 1, j, k)] / hx / hx +
           data[st][index(i, j - 1, k)] - 2 * data[st][index(i, j, k)] + data[st][index(i, j + 1, k)] / hy / hy +
           data[st][index(i, j, k - 1)] - 2 * data[st][index(i, j, k)] + data[st][index(i, j, k + 1)] / hz / hz;
}


void assign_boundaries(double **data, int stage) {
    // x
    for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nz; k++) {
            data[stage][index(0, j, k)] = data[stage][index(Nx - 2, j, k)];
            data[stage][index(Nx - 1, j, k)] = data[stage][index(1, j, k)];
        }
    }

    // y
    for (int i = 0; i < Nx; i++) {
        for (int k = 0; k < Nz; k++) {
            data[stage][index(i, 0, k)] = 0;
            data[stage][index(i, Ny - 1, k)] = 0;
        }
    }

    // z
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            data[stage][index(i, j, 0)] = data[stage][index(i, j, Nz - 2)];
            data[stage][index(i, j, Nz - 1)] = data[stage][index(i, j, 1)];
        }
    }
}

void init(double **data) {

    // t = t0
    for (int idx = 0; idx < Nx * Ny * Nz; idx++) {
        auto [i, j, k] = idx_to_ijk(idx);
        data[0][idx] = func_phi(real_x(i), real_y(j), real_z(k));
    }

    // t = t1
    for (int idx = 0; idx < Nx * Ny * Nz; idx++) {
        auto [i, j, k] = idx_to_ijk(idx);

        if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1 || k == 0 || k == Nz - 1)
            continue;

        data[1][idx] = data[0][idx] + tau * tau / 2 * laplace((const double**)data, 0, i, j, k);
    }

    assign_boundaries(data, 1);
}

void calculate(double **data) {
    for (int t = 2; t < K; t++) {
        for (int idx = 0; idx < Nx * Ny * Nz; idx++) {
            auto [i, j, k] = idx_to_ijk(idx);

            if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1 || k == 0 || k == Nz - 1)
                continue;

            data[t % num_stages][idx] = 2 * data[(t - 1) % num_stages][idx] - data[(t - 2) % num_stages][idx]
                                        + tau * tau * a_2 * laplace((const double**)data, (t - 1) % num_stages, i, j, k);
            assign_boundaries(data, t % num_stages);
        }
        std::cout << "t = " << t << std::endl;
    }
}


int main(int argc, char **argv) {

    N = (argc > 1) ? std::atoi(argv[1]) : 128;
    Nx = Ny = Nz = N;

    // Lx = Ly = Lz = (argv[2] == "pi") ? M_PI : std::atoi(argv[2]);
    Lx = Ly = Lz = 1;

    hx = Lx / (Nx - 1);
    hy = Ly / (Ny - 1);
    hz = Lz / (Nz - 1);
    
    Nx += 1;
    Nz += 1;

    // std::cout << N << " " << Nx << " " << Ny << " " << Nz << std::endl;

    double *data[num_stages];
    for (int st = 0; st < num_stages; st++)
        data[st] = new double[Nx * Ny * Nz];

    init(data);
    calculate(data);


    // for (int n = 0; n < 2; n++) {
    //     for (int j = 0; j < Ny; j++) {
    //         for (int i = 0; i < Nx; i++) {
    //             int idx = index(i, j, 0);
    //             std::cout << data[n][idx] << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << std::vector<double>(data[0], data[0] + Nx * Ny * Nz).size() << std::endl;
    // std::cout << data[0][10] << std::endl;

    // auto [i, j, k] = idx_to_ijk(20000);


    return 0;
}