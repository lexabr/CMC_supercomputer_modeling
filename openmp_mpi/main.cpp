#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <omp.h>
#include <mpi.h>


int N, Nx, Ny, Nz;
double Lx, Ly, Lz;
double hx, hy, hz;

const double a_2 = 1/9;
const double tau = 0.001;
const int K = 20;
const int num_stages = 3;
int ndim = 3;

std::string prog_postfix;

int rank, nproc;
int i_min, j_min, k_min, i_max, j_max, k_max;


inline double func_u(double x, double y, double z, double t) {
    double a_t = M_PI / 3 * sqrt(4 / pow(Lx, 2) + 1 / pow(Ly, 2) + 4 / pow(Lz, 2));

    return sin(2 * M_PI * x / Lx) * sin(M_PI * y / Ly + M_PI) * sin(2 * M_PI * z / Lz + 2 * M_PI) * cos(a_t * t + M_PI);
}

inline double func_phi(double x, double y, double z) {
    return func_u(x, y, z, 0);
}

inline double real_x(int i) {return i * hx;}
inline double real_y(int j) {return j * hy;}
inline double real_z(int k) {return k * hz;}
inline double real_t(int t) {return t * tau;}



// utils
inline int index(int i, int j, int k) { // возвращает индекс элемента внутри блока
    return Ny * Nx * k + Nx * j + i;
}

double laplace(const double **data, int st, int i, int j, int k) {

    return (data[st][index(i - 1, j, k)] - 2 * data[st][index(i, j, k)] + data[st][index(i + 1, j, k)]) / (hx * hx) +
           (data[st][index(i, j - 1, k)] - 2 * data[st][index(i, j, k)] + data[st][index(i, j + 1, k)]) / (hy * hy) +
           (data[st][index(i, j, k - 1)] - 2 * data[st][index(i, j, k)] + data[st][index(i, j, k + 1)]) / (hz * hz);
}

// double* real_u() {
//     double *u_data = new double[Nx * Ny * Nz];

//     #pragma omp parallel for collapse(3)
//     for (int i = 0; i < Nx; i++)
//         for (int j = 0; j < Ny; j++)
//             for (int k = 0; k < Nz; k++)
//                 u_data[index(i, j, k)] = func_u(real_x(i), real_y(j), real_z(k), real_t(K - 1));
    
//     return u_data;
// }


// class Error {
// private:
//     double rmse;
//     double max_abs;
//     double *grid_errs;

// public:
//     Error() {
//         rmse = 0;
//         max_abs = 0;
//         grid_errs = new double[Nx * Ny * Nz];
//     }

//     ~Error() {
//         delete[] grid_errs;
//     }

//     double get_rmse()       const {return this->rmse;}
//     double get_max_abs()    const {return this->max_abs;}
//     double* get_grid_errs() const {return this->grid_errs;}

//     void calc_stage_errs(const double **data, int st) {

//         double mse = 0;
//         double max_abs = 0;
        
//         #pragma omp parallel for collapse(3) reduction(+:mse) reduction(max: max_abs)
//         for (int i = 0; i < Nx; i++)
//             for (int j = 0; j < Ny; j++)
//                 for (int k = 0; k < Nz; k++) {
//                     double diff = func_u(real_x(i), real_y(j), real_z(k), real_t(st)) - data[st % num_stages][index(i, j, k)];
//                     mse += pow(diff, 2);

//                     max_abs = std::max(max_abs, fabs(diff));
//                 }


//         mse /= Nx * Ny * Nz;

//         this->rmse = sqrt(mse);
//         this->max_abs = max_abs;

//         this->write_stage_errs(st, st == 0);
//     }

//     void write_stage_errs(int st, bool truncate = false) {
//         std::ofstream file;
//         std::string filename = "./results/stage_errors" + prog_postfix + ".csv";

//         if (truncate) {
//             file.open(filename);
//             file << "t,RMSE,max_abs" << std::endl;
//         } else {
//             file.open(filename, std::ios_base::app);
//         }
        
//         file << st << "," << this->rmse << "," << this->max_abs << std::endl;
//         file.close();
//     }

//     void calc_grid_errs(const double **data, int st) {
//         #pragma omp parallel for collapse(3)
//         for (int i = 0; i < Nx; i++)
//             for (int j = 0; j < Ny; j++)
//                 for (int k = 0; k < Nz; k++)
//                     this->grid_errs[index(i, j, k)] = func_u(real_x(i), real_y(j), real_z(k), real_t(st)) - data[st % num_stages][index(i, j, k)];
//     }
// };


class Timer {
public:
    double start_time;
    double end_time;
    double overall_dur;
    double calc_dur;
    double calc_errors_time;

    void start() {this->start_time = MPI_Wtime();}
    void end()   {
        this->end_time = MPI_Wtime();
        this->calc();
    }
    void add_err_time(double et) {this->calc_errors_time += et;}

    void overall_duration() {this->overall_dur = this->end_time - this->start_time;}
    void calculation_duration() {this->calc_dur = this->end_time - this->start_time - this->calc_errors_time;}
    void calc() {
        this->overall_duration();
        this->calculation_duration();
    }

    double get_overall_duration() const {return this->overall_dur;}
    double get_calculation_duration() const {return this->calc_dur;}
    double get_error_calculation_duration() const {return this->calc_errors_time;}

    void write_duration() {
        std::ofstream file;
        file.open("./results/info" + prog_postfix + ".txt", std::ios_base::app);

        file << "N = " << N << std::endl;
        file << "L = " << Lx << std::endl;
        file << "K = " << K << std::endl;
        file << "tau = " << tau << std::endl;
        file << "Total dur (s): " << this->get_overall_duration() << std::endl;
        file << "Calculations dur (s): " << this->get_calculation_duration() << std::endl;
        file << "Error calculations & file writing dur (s): " << this->get_error_calculation_duration() << std::endl << std::endl;
    }
};


// void save_grid_values(const double *data, const std::string filename) {
//     std::ofstream file;
//     file.open(filename);

//     for (int k = 0; k < Nz - 1; k++)
//         for (int j = 0; j < Ny; j++) {
//             for (int i = 0; i < Nx - 1; i++) {
//                 file << data[index(i, j, k)];
//                 if (i < Nx - 2)
//                     file << ",";
//             }
//             file << std::endl;
//         }
//     file.close();
// }


void send_recv_right(double *data, int axis, int bound_cond, MPI_Comm& comm_cart, int rank_prev, int rank_next, bool is_first, bool is_last) {
    int size, dim1, dim2;

    switch (axis) {
        case 0:
            dim1 = Ny;
            dim2 = Nz;
            break;
        case 1:
            dim1 = Nx;
            dim2 = Nz;
            break;
        case 2:
            dim1 = Nx;
            dim2 = Ny;    
        default:
            break;
    }
    size = dim1 * dim2;

    if (is_first && is_last) // equal to dims[axis] == 1
        return;

    double send_buf[size], recv_buf[size];
    MPI_Status comm_status;

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++) {
            int idx;
            switch (axis) {
                case 0:
                    idx = index((is_last && bound_cond == 2) ? Nx - 3 : Nx - 2, i, j);
                    break;
                case 1:
                    idx = index(i, (is_last && bound_cond == 2) ? Ny - 3 : Ny - 2, j);
                    break;
                case 2:
                    idx = index(i, j, (is_last && bound_cond == 2) ? Nz - 3 : Nz - 2);
                    break;
                default:
                    break;
            }
            send_buf[dim1*j+i] = data[idx];
        }

    MPI_Sendrecv(send_buf, size, MPI_DOUBLE, rank_next, 1,
                 recv_buf, size, MPI_DOUBLE, rank_prev, 1,
                 comm_cart, &comm_status);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++) {
            int idx;
            switch (axis) {
                case 0:
                    idx = index(0, i, j);
                    break;
                case 1:
                    idx = index(i, 0, j);
                    break;
                case 2:
                    idx = index(i, j, 0);
                    break;
                default:
                    break;
            }
            data[idx] = recv_buf[dim1*j+i];
        }
}


void send_recv_left(double *data, int axis, int bound_cond, MPI_Comm& comm_cart, int rank_prev, int rank_next, bool is_first, bool is_last) {
    int size, dim1, dim2;

    switch (axis) {
        case 0:
            dim1 = Ny;
            dim2 = Nz;
            break;
        case 1:
            dim1 = Nx;
            dim2 = Nz;
            break;
        case 2:
            dim1 = Nx;
            dim2 = Ny;    
        default:
            break;
    }
    size = dim1 * dim2;

    if (is_first && is_last) // equal to dims[axis] == 1
        return;

    double send_buf[size], recv_buf[size];
    MPI_Status comm_status;

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++) {
            int idx;
            switch (axis) {
                case 0:
                    idx = index((is_first && bound_cond == 2) ? 2 : 1, i, j);
                    break;
                case 1:
                    idx = index(i, (is_first && bound_cond == 2) ? 2 : 1, j);
                    break;
                case 2:
                    idx = index(i, j, (is_first && bound_cond == 2) ? 2 : 1);
                    break;
                default:
                    break;
            }
            send_buf[dim1*j+i] = data[idx];
        }

    MPI_Sendrecv(send_buf, size, MPI_DOUBLE, rank_prev, 1,
                 recv_buf, size, MPI_DOUBLE, rank_next, 1,
                 comm_cart, &comm_status);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++) {
            int idx;
            switch (axis) {
                case 0:
                    idx = index(Nx - 1, i, j);
                    break;
                case 1:
                    idx = index(i, Ny - 1, j);
                    break;
                case 2:
                    idx = index(i, j, Nz - 1);
                    break;
                default:
                    break;
            }
            data[idx] = recv_buf[dim1*j+i];
        }
}

void send_recv(double *data, MPI_Comm& comm_cart, int bound_cond[], int rank_prev[], int rank_next[], bool is_first[], bool is_last[]) {

    #pragma omp parallel for
    for (int d = 0; d < ndim; d++) {
        send_recv_right(data, d, bound_cond[d], comm_cart, rank_prev[d], rank_next[d], is_first[d], is_last[d]);
        send_recv_left(data, d, bound_cond[d], comm_cart, rank_prev[d], rank_next[d], is_first[d], is_last[d]);
    }
}


void assign_boundaries(double **data, int stage, bool is_first[], bool is_last[]) {
    // x
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nz; k++) {
            if (is_first[0] && is_last[0]) {
                data[stage][index(1, j, k)] = data[stage][index(Nx - 3, j, k)];
                data[stage][index(Nx - 2, j, k)] = data[stage][index(2, j, k)];
            } else if (is_first[0]) {
                data[stage][index(1, j, k)] = data[stage][index(0, j, k)]; // с учетом обменных областей
            } else if (is_last[0]) {
                data[stage][index(Nx - 2, j, k)] = data[stage][index(Nx - 1, j, k)]; // с учетом обменных областей
            }
        }
    }

    // y
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < Nx; i++) {
        for (int k = 0; k < Nz; k++) {
            if (is_first[1])
                data[stage][index(i, 1, k)] = 0;
            else if (is_last[1])
                data[stage][index(i, Ny - 2, k)] = 0;
        }
    }

    // z
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            if (is_first[2] && is_last[2]) {
                data[stage][index(i, j, 1)] = data[stage][index(i, j, Nz - 3)];
                data[stage][index(i, j, Nz - 2)] = data[stage][index(i, j, 2)];
            } else if (is_first[2]) {
                data[stage][index(i, j, 1)] = data[stage][index(i, j, 0)]; // с учетом обменных областей
            } else if (is_last[2]) {
                data[stage][index(i, j, Nz - 2)] = data[stage][index(i, j, Nz - 1)]; // с учетом обменных областей
            }
        }
    }
}

void init(double **data, MPI_Comm& comm_cart, int bound_cond[], int rank_prev[], int rank_next[],
            bool is_first[], bool is_last[], Timer& timer) {//, Error &err, Timer &timer) {

    // t = t0
    #pragma omp parallel for collapse(3)
    for (int i = 1; i < Nx - 1; i++)
        for (int j = 1; j < Ny - 1; j++)
            for (int k = 1; k < Nz - 1; k++)
                data[0][index(i, j, k)] = func_phi(real_x(i_min + i - 1), real_y(j_min + j - 1), real_z(k_min + k - 1)); // -1 из-за левой обменной области

    double err_st = MPI_Wtime();
    // err.calc_stage_errs((const double **)data, 0);
    timer.add_err_time(MPI_Wtime() - err_st);

    // t = t1
    #pragma omp parallel for collapse(3)
    for (int i = 1; i < Nx - 1; i++)
        for (int j = 1; j < Ny - 1; j++)
            for (int k = 1; k < Nz - 1; k++) {
                if ((i == 1 && is_first[0]) || (i == Nx - 2 && is_last[0]))
                    continue;
                if ((j == 1 && is_first[1]) || (j == Ny - 2 && is_last[1]))
                    continue;
                if ((k == 1 && is_first[2]) || (k == Nz - 2 && is_last[2]))
                    continue;

                data[1][index(i, j, k)] = data[0][index(i, j, k)] + tau * tau * a_2 / 2 * laplace((const double**)data, 0, i, j, k);
            }

    send_recv(data[1], comm_cart, bound_cond, rank_prev, rank_next, is_first, is_last);

    assign_boundaries(data, 1, is_first, is_last);

    err_st = MPI_Wtime();
    // err.calc_stage_errs((const double **)data, 1);
    timer.add_err_time(MPI_Wtime() - err_st);
}

void calculate(double **data, MPI_Comm& comm_cart, int bound_cond[], int rank_prev[], int rank_next[],
                bool is_first[], bool is_last[], Timer& timer) {
    // Error err;

    init(data, comm_cart, bound_cond, rank_prev, rank_next, is_first, is_last, timer);


    for (int t = 2; t < K; t++) {
        #pragma omp parallel for collapse(3)
        for (int i = 1; i < Nx - 1; i++)
            for (int j = 1; j < Ny - 1; j++)
                for (int k = 1; k < Nz - 1; k++) {
                    if ((i == 1 && is_first[0]) || (i == Nx - 2 && is_last[0]))
                        continue;
                    if ((j == 1 && is_first[1]) || (j == Ny - 2 && is_last[1]))
                        continue;
                    if ((k == 1 && is_first[2]) || (k == Nz - 2 && is_last[2]))
                        continue;

                    data[t % num_stages][index(i, j, k)] = 2 * data[(t - 1) % num_stages][index(i, j, k)] - data[(t - 2) % num_stages][index(i, j, k)]
                                                            + tau * tau * a_2 * laplace((const double**)data, (t - 1) % num_stages, i, j, k);
                }
        
        assign_boundaries(data, t % num_stages, is_first, is_last);

        double err_st = MPI_Wtime();
        // err.calc_stage_errs((const double**)data, t);
        timer.add_err_time(MPI_Wtime() - err_st);
    }

    // double err_st = omp_get_wtime();
    // err.calc_grid_errs((const double**)data, K - 1);
    // if (write_last_stage) {
    //     save_grid_values((const double*)err.get_grid_errs(), "./results/grid_errs" + prog_postfix + ".csv");
    //     save_grid_values((const double*)data[(K - 1) % num_stages], "./results/u_data" + prog_postfix + ".csv");

    //     double *real_u_data = real_u();
    //     save_grid_values((const double*)real_u_data, "./results/u_real" + prog_postfix + ".csv");
    //     delete[] real_u_data;
    // }
    
    // timer.add_err_time(omp_get_wtime() - err_st);
}



int main(int argc, char **argv) {
    N = (argc > 1) ? std::atoi(argv[1]) : 128;
    Lx = Ly = Lz = ((argc > 2) && (std::atoi(argv[2]) == 2)) ? M_PI : 1;
    bool need_write = (argc > 3) ? (bool)std::atoi(argv[3]) : false;
    int num_threads = (argc > 4) ? std::atoi(argv[4]) : 1;

    hx = Lx / (N - 1);
    hy = Ly / (N - 1);
    hz = Lz / (N - 1);


    int boundary_conditions[ndim];
    boundary_conditions[0] = 2;
    boundary_conditions[1] = 1;
    boundary_conditions[2] = 2;

    int dims[ndim];
    dims[0] = 0;
    dims[1] = 0;
    dims[2] = 0;

    int periods[ndim];
    for (int i = 0; i < ndim; i++)
        periods[i] = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    Timer timer_proc, timer_all;
    timer_proc.start();

    // вычисление оптимальных размерностей решетки
    MPI_Dims_create(nproc, ndim, dims);

    // if (rank == 0) {
    //     for (int i = 0; i < ndim; i++)
    //         std::cout << dims[i];
    //     std::cout << std::endl << std::endl;;
    // }

    // число узлов сетки внутри блока
    int nodes[ndim];
    for (int i = 0; i < ndim; i++) {
        nodes[i] = N / dims[i];
    }

    // создание топологии
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, 0, &comm_cart);

    // координаты блока внутри решетки
    int coords[ndim];
    MPI_Cart_coords(comm_cart, rank, ndim, coords);

    // ранги соседей по каждой оси для обменов
    int rank_prev[ndim], rank_next[ndim];
    for (int i = 0; i < ndim; i++)
        MPI_Cart_shift(comm_cart, i, 1, &rank_prev[i], &rank_next[i]);

    bool is_first[ndim], is_last[ndim];
    for (int i = 0; i < ndim; i++) {
        is_first[i] = coords[i] == 0;
        is_last[i] = coords[i] == dims[i] - 1;
    }

    // минимальные и максимальные индексы внутри общей сетки
    i_min = coords[0] * nodes[0];
    j_min = coords[1] * nodes[1];
    k_min = coords[2] * nodes[2];
    i_max = (coords[0] + 1) * nodes[0] - ((coords[0] == dims[0] - 1) && (boundary_conditions[0] == 2) ? 0 : 1); // коррекция на случай периодических граничных условий
    j_max = (coords[1] + 1) * nodes[1] - ((coords[1] == dims[1] - 1) && (boundary_conditions[1] == 2) ? 0 : 1); // коррекция на случай периодических граничных условий
    k_max = (coords[2] + 1) * nodes[2] - ((coords[2] == dims[2] - 1) && (boundary_conditions[2] == 2) ? 0 : 1); // коррекция на случай периодических граничных условий

    // объём выделяемой памяти (+3 = +1 для точек функции и +2 для обменных областей слева и справа)
    Nx = i_max - i_min + 3;
    Ny = j_max - j_min + 3;
    Nz = k_max - k_min + 3;

    prog_postfix = "_N_" + std::to_string(N) + "_L_" + ((Lx == 1) ? std::to_string(1) : "Pi") + "_P_" + std::to_string(nproc) + "_Th_" + std::to_string(num_threads);

    double *data[num_stages];
    for (int st = 0; st < num_stages; st++)
        data[st] = new double[Nx * Ny * Nz];
    // std::cout << rank << " " << Nx * Ny * Nz << std::endl;

    // std::cout << "Process: " << rank << std::endl;
    // for (int j = 0; j < Ny; j++) {
    //     for (int i = 0; i < Nx; i++) {
    //         std::cout << data[1][index(i, j, 1)] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
    calculate(data, comm_cart, boundary_conditions, rank_prev, rank_next, is_first, is_last, timer_proc);

    timer_proc.end();

    timer_all.overall_dur = 0;
    timer_all.calc_dur = 0;
    timer_all.calc_errors_time = 0;
    MPI_Reduce(&timer_proc.overall_dur, &timer_all.overall_dur, 1, MPI_DOUBLE, MPI_MAX, 0, comm_cart);
    MPI_Reduce(&timer_proc.calc_dur, &timer_all.calc_dur, 1, MPI_DOUBLE, MPI_MAX, 0, comm_cart);
    MPI_Reduce(&timer_proc.calc_errors_time, &timer_all.calc_errors_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm_cart);

    // std::cout << timer_proc.get_overall_duration() << " " << timer_proc.get_calculation_duration() << " " << timer_proc.get_error_calculation_duration() << std::endl;
    if (rank == 0)
        timer_all.write_duration();
    //     std::cout << "gg " << timer_all.get_overall_duration() << " " << timer_all.get_calculation_duration() << " " << timer_all.get_error_calculation_duration() << std::endl;

    MPI_Finalize();
    
    // timer.write_duration();

    return 0;
}