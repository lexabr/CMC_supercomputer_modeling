#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <ctime>

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

double* real_u() {
    double *u_data = new double[Nx * Ny * Nz];

    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
                u_data[index(i, j, k)] = func_u(real_x(i), real_y(j), real_z(k), real_t(K - 1));
    
    return u_data;
}


class Error {
private:
    double rmse;
    double max_abs;
    double *grid_errs;

public:
    Error() {
        rmse = 0;
        max_abs = 0;
        grid_errs = new double[Nx * Ny * Nz];
    }

    ~Error() {
        delete[] grid_errs;
    }

    double get_rmse()       const {return this->rmse;}
    double get_max_abs()    const {return this->max_abs;}
    double* get_grid_errs() const {return this->grid_errs;}

    void calc_stage_errs(const double **data, int st) {

        double mse = 0;
        double max_abs = 0;

        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++) {
                    double diff = func_u(real_x(i), real_y(j), real_z(k), real_t(st)) - data[st % num_stages][index(i, j, k)];
                    mse += pow(diff, 2);

                    double abs_diff = fabs(diff);
                    if (abs_diff > max_abs)
                        max_abs = abs_diff;
                }


        mse /= Nx * Ny * Nz;

        this->rmse = sqrt(mse);
        this->max_abs = max_abs;

        this->write_stage_errs(st, st == 0);
    }

    void write_stage_errs(int st, bool truncate = false) {
        std::ofstream file;
        std::string filename = "./results/stage_errors_N_" + std::to_string(N) + "_.csv";

        if (truncate) {
            file.open(filename);
            file << "t,RMSE,max_abs" << std::endl;
        } else {
            file.open(filename, std::ios_base::app);
        }
        
        file << st << "," << this->rmse << "," << this->max_abs << std::endl;
        file.close();
    }

    void calc_grid_errs(const double **data, int st) {
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++)
                    this->grid_errs[index(i, j, k)] = func_u(real_x(i), real_y(j), real_z(k), real_t(st)) - data[st % num_stages][index(i, j, k)];
    }
};


class Timer {
public:
    std::clock_t start_time;
    std::clock_t end_time;
    std::clock_t calc_errors_time = 0;

public:

    void start() {this->start_time = std::clock();}
    void end()   {this->end_time = std::clock();}
    void add_err_time(std::clock_t et) {this->calc_errors_time += et;}

    double overall_duration()           const {return (this->end_time - this->start_time) / (double)CLOCKS_PER_SEC;}
    double calculation_duration()       const {return (this->end_time - this->start_time - this->calc_errors_time) / (double)CLOCKS_PER_SEC;}
    double error_calculation_duration() const {return this->calc_errors_time / (double)CLOCKS_PER_SEC;}

    void write_duration() {
        std::ofstream file;
        file.open("./results/info_N_" + std::to_string(N) + "_.txt");

        file << "N = " << N << std::endl;
        file << "L = " << Lx << std::endl;
        file << "K = " << K << std::endl;
        file << "tau = " << tau << std::endl;
        file << "Total dur (s): " << this->overall_duration() << std::endl;
        file << "Calculations dur (s): " << this->calculation_duration() << std::endl;
        file << "Error calculations & file writing dur (s): " << this->error_calculation_duration() << std::endl;
    }
};


void save_grid_values(const double *data, const std::string filename) {
    std::ofstream file;
    file.open(filename);

    for (int k = 0; k < Nz - 1; k++)
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx - 1; i++) {
                file << data[index(i, j, k)];
                if (i < Nx - 2)
                    file << ",";
            }
            file << std::endl;
        }
    file.close();
}


double laplace(const double **data, int st, int i, int j, int k) {

    return (data[st][index(i - 1, j, k)] - 2 * data[st][index(i, j, k)] + data[st][index(i + 1, j, k)]) / (hx * hx) +
           (data[st][index(i, j - 1, k)] - 2 * data[st][index(i, j, k)] + data[st][index(i, j + 1, k)]) / (hy * hy) +
           (data[st][index(i, j, k - 1)] - 2 * data[st][index(i, j, k)] + data[st][index(i, j, k + 1)]) / (hz * hz);
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

void init(double **data, Error &err, Timer &timer) {

    // t = t0
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
                data[0][index(i, j, k)] = func_phi(real_x(i), real_y(j), real_z(k));

    clock_t err_st = std::clock();
    err.calc_stage_errs((const double **)data, 0);
    timer.add_err_time(std::clock() - err_st);

    // t = t1
    for (int i = 1; i < Nx - 1; i++)
        for (int j = 1; j < Ny - 1; j++)
            for (int k = 1; k < Nz - 1; k++)
                data[1][index(i, j, k)] = data[0][index(i, j, k)] + tau * tau * a_2 / 2 * laplace((const double**)data, 0, i, j, k);

    assign_boundaries(data, 1);

    err_st = std::clock();
    err.calc_stage_errs((const double **)data, 1);
    timer.add_err_time(std::clock() - err_st);
}

void calculate(double **data, Timer &timer, bool write_last_stage = false) {
    Error err;

    init(data, err, timer);

    for (int t = 2; t < K; t++) {
        for (int i = 1; i < Nx - 1; i++)
            for (int j = 1; j < Ny - 1; j++)
                for (int k = 1; k < Nz - 1; k++) 
                    data[t % num_stages][index(i, j, k)] = 2 * data[(t - 1) % num_stages][index(i, j, k)] - data[(t - 2) % num_stages][index(i, j, k)]
                                                            + tau * tau * a_2 * laplace((const double**)data, (t - 1) % num_stages, i, j, k);
        
        assign_boundaries(data, t % num_stages);

        clock_t err_st = std::clock();
        err.calc_stage_errs((const double**)data, t);
        timer.add_err_time(std::clock() - err_st);
    }

    clock_t err_st = std::clock();
    err.calc_grid_errs((const double**)data, K - 1);
    if (write_last_stage) {
        save_grid_values((const double*)err.get_grid_errs(), "./results/grid_errs_N_" + std::to_string(N) + "_.csv");
        save_grid_values((const double*)data[(K - 1) % num_stages], "./results/u_data_N_" + std::to_string(N) + "_.csv");

        double *real_u_data = real_u();
        save_grid_values((const double*)real_u_data, "./results/u_real_N_" + std::to_string(N) + "_.csv");
        delete[] real_u_data;
    }
    
    timer.add_err_time(std::clock() - err_st);
}


int main(int argc, char **argv) {
    Timer timer;
    timer.start();

    N = (argc > 1) ? std::atoi(argv[1]) : 128;
    Nx = Ny = Nz = N;

    // Lx = Ly = Lz = (argv[2] == "pi") ? M_PI : std::atoi(argv[2]);
    Lx = Ly = Lz = 1;

    hx = Lx / (N - 1);
    hy = Ly / (N - 1);
    hz = Lz / (N - 1);
    
    Nx += 1;
    Nz += 1;


    double *data[num_stages];
    for (int st = 0; st < num_stages; st++)
        data[st] = new double[Nx * Ny * Nz];

    calculate(data, timer, false);
    
    timer.end();
    timer.write_duration();

    return 0;
}