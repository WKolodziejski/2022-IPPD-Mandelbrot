#include <iostream>
#include <complex>
#include<string>

using std::string;
using std::cout;

constexpr auto max_row = 22, max_column = 78, max_iteration = 20;

int main() {
    char output[max_row][max_column];

    #pragma omp parallel default(none) shared(output, cout)
    {
        #pragma omp for
        for (auto row = 0; row < max_row; ++row) {

            #pragma omp parallel default(none) shared(output, cout, row)
            {
                #pragma omp for
                for (auto column = 0; column < max_column; ++column) {
                    std::complex<float> z, c = {
                            (float) column * 2.0f / max_column - 1.5f,
                            (float) row * 2.0f / max_row - 1
                    };

                    int i = 0;
                    while (abs(z) < 2 && ++i < max_iteration)
                        z = pow(z, 2) + c;

                    output[row][column] = (i == max_iteration ? '#' : '.');
                    cout << (i == max_iteration ? '#' : '.');
                }
            }
            cout << '\n';
        }
    }

    for (auto & row : output) {
        for (char column : row) {
            std::cout << column;
        }
        std::cout << '\n';
    }

    return 0;
}
