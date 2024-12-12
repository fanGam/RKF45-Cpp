#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <GL/glut.h>

class RKF45 {
private:
    // Коэффициенты метода Рунге-Кутты-Фельберга
    const double a2 = 1.0 / 4.0;
    const double a3 = 3.0 / 8.0;
    const double a4 = 12.0 / 13.0;
    const double a5 = 1.0;
    const double a6 = 1.0 / 2.0;
    
    const double b21 = 1.0 / 4.0;
    const double b31 = 3.0 / 32.0;
    const double b32 = 9.0 / 32.0;
    const double b41 = 1932.0 / 2197.0;
    const double b42 = -7200.0 / 2197.0;
    const double b43 = 7296.0 / 2197.0;
    const double b51 = 439.0 / 216.0;
    const double b52 = -8.0;
    const double b53 = 3680.0 / 513.0;
    const double b54 = -845.0 / 4104.0;
    const double b61 = -8.0 / 27.0;
    const double b62 = 2.0;
    const double b63 = -3544.0 / 2565.0;
    const double b64 = 1859.0 / 4104.0;
    const double b65 = -11.0 / 40.0;

    const double c1 = 16.0 / 135.0;
    const double c3 = 6656.0 / 12825.0;
    const double c4 = 28561.0 / 56430.0;
    const double c5 = -9.0 / 50.0;
    const double c6 = 2.0 / 55.0;

    const double d1 = 25.0 / 216.0;
    const double d3 = 1408.0 / 2565.0;
    const double d4 = 2197.0 / 4104.0;
    const double d5 = -1.0 / 5.0;

    // Функция для вычисления производной
    // Химическая реакция
    double him(double t, double y) {
        double k = 0.05; // Скорость химической реакции
        return -k * y * y;
    }

    // Падение тела
    double fall(double t, double y) {
        double g = 9.8; // Гравитационная постоянная
        double k = 0.1; // Коэфф сопротивления воздуха
        double m = 1; // Масса тела
        return g - (k / m) * y * y;
    }

    // Модель роста популяции
    double popul(double t, double y) {
        return 0.1 * y * (1 - y / 100.);
    }

    // Радиоактивный распад
    double radio(double t, double quantity) {
        double half_life = 5730;
        return -log(2) / half_life * quantity;
    }

public:
    // Метод решения ОДУ с адаптивным шагом
    std::vector<std::pair<double, double>> solve(
        double t0,     // начальное время
        double y0,     // начальное значение
        double tf,     // конечное время
        double h0,     // начальный шаг
        double tol,    // допустимая погрешность,
        int met        // метод решения
    ) {
        std::vector<std::pair<double, double>> solution;
        solution.push_back({ t0, y0 });

        double t = t0;
        double y = y0;
        double h = h0;

        while (t < tf) {
            // Ограничение максимального шага
            if (t + h > tf) {
                h = tf - t;
            }

            // Вычисление коэффициентов Рунге-Кутты
            double k1;
            double k2;
            double k3;
            double k4;
            double k5;
            double k6;
            if (met == 0) {
                k1 = h * radio(t, y);
                k2 = h * radio(t + a2 * h, y + b21 * k1);
                k3 = h * radio(t + a3 * h, y + b31 * k1 + b32 * k2);
                k4 = h * radio(t + a4 * h, y + b41 * k1 + b42 * k2 + b43 * k3);
                k5 = h * radio(t + a5 * h, y + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4);
                k6 = h * radio(t + a6 * h, y + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5);
            }
            if (met == 1) {
                k1 = h * popul(t, y);
                k2 = h * popul(t + a2 * h, y + b21 * k1);
                k3 = h * popul(t + a3 * h, y + b31 * k1 + b32 * k2);
                k4 = h * popul(t + a4 * h, y + b41 * k1 + b42 * k2 + b43 * k3);
                k5 = h * popul(t + a5 * h, y + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4);
                k6 = h * popul(t + a6 * h, y + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5);
            }
            if (met == 2) {
                k1 = h * fall(t, y);
                k2 = h * fall(t + a2 * h, y + b21 * k1);
                k3 = h * fall(t + a3 * h, y + b31 * k1 + b32 * k2);
                k4 = h * fall(t + a4 * h, y + b41 * k1 + b42 * k2 + b43 * k3);
                k5 = h * fall(t + a5 * h, y + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4);
                k6 = h * fall(t + a6 * h, y + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5);
            }
            if (met == 3) {
                k1 = h * him(t, y);
                k2 = h * him(t + a2 * h, y + b21 * k1);
                k3 = h * him(t + a3 * h, y + b31 * k1 + b32 * k2);
                k4 = h * him(t + a4 * h, y + b41 * k1 + b42 * k2 + b43 * k3);
                k5 = h * him(t + a5 * h, y + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4);
                k6 = h * him(t + a6 * h, y + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5);
            }

            // Вычисление значений по 4 и 5 порядку точности
            double y4 = y + d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5;
            double y5 = y + c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6;

            // Оценка локальной погрешности
            double error = std::abs(y5 - y4);
            double scaled_error = error / h;

            // Адаптация шага
            if (scaled_error <= tol) {
                t += h;
                y = y5;
                solution.push_back({ t, y });
            }

            // Пересчет шага
            double h_new = 0.9 * h * std::pow(tol / scaled_error, 0.2);
            h = std::min(std::max(h_new, 0.1 * h), 5.0 * h);
        }

        return solution;
    }
};

int main() {
    RKF45 solver;
    // с начальными условиями: y(0), t(1), t(1), step, sq
    //auto solution = solver.solve(t(0), y(0), t(1), step, eq);
    // Радиоактивный распад
    //auto solution = solver.solve(0.0, 100.0, 20000.0, 100, 1e-6, 0);
    // Модель роста популяции
    //auto solution = solver.solve(0.0, 10.0, 100.0, 0.1, 1e-6, 1);
    // Падение тела
    //auto solution = solver.solve(0.0, 0.0, 10.0, 0.1, 1e-6, 2);
    // Химическая реакция
    auto solution = solver.solve(0.0, 100.0, 20.0, 0.1, 1e-6, 3);
    
    // Вывод результатов
    for (const auto& point : solution) {
        //std::cout << "t = " << point.first
            //<< ", y = " << point.second << std::endl;
        std::cout << point.first
            << " " << point.second << std::endl;
    }

    return 0;
}