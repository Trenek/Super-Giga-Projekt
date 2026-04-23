#include <capd/capdlib.h>

#define gamma 2
#define beta 4
#define k 1

typedef capd::autodiff::Node num;

template <size_t N>
static num f(num x0, num xN, num n) {
    return - gamma * x0 + beta * (xN ^ k) / (1 + exp(n * log(xN)));
}

template <size_t N>
static double x(double j) {
    return cos(j * std::numbers::pi / N);
}

template <size_t N>
static double c(double i) {
    return 1 + (i == 0) + (i == N);
}

template <size_t N>
static std::array<std::array<double, N + 1>, N + 1> generateChebyshev() {
    std::array<std::array<double, N + 1>, N + 1> result;

    for (int i = 0; i < N + 1; i += 1) {
        for (int j = 0; j < N + 1; j += 1) {
            result[i][j] =
                (i == 0 && j == 0) ?   (2.0 * N * N + 1) / 6 :
                (i == N && j == N) ? - (2.0 * N * N + 1) / 6 :
                (i == j)           ? - x<N>(i) / (2 * (1 - pow(x<N>(i), 2))) :
                                       c<N>(i) * pow(-1, i + j) / (c<N>(j) * (x<N>(i) - x<N>(j)));
            result[i][j] *= 2;
        }
    }

    return result;
}

template <size_t N>
static void applyChebyshev(num in[], num out[]) {
    static std::array<std::array<double, N + 1>, N + 1> M = generateChebyshev<N>();

    for (size_t i = 1; i <= N; i += 1) {
        out[i] = M[i][0] * in[0];
        for (size_t j = 1; j <= N; j += 1) {
            out[i] += M[i][j] * in[j];
        }
    }
}

template <size_t N>
void mackeyGlass(num &time, num in[], int dimIn, num out[], int dimOut, num param[], int noParam) {
    out[0] = f<N>(in[0], in[dimIn - 1], param[0]);
    applyChebyshev<N>(in, out);
}

#undef gamma
#undef beta
#undef k
