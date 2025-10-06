#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tinyexpr/tinyexpr.h"

char funcao_str[256];

double funcao(double x) {
    te_variable vars[] = { {"x", &x} };
    int err;
    te_expr *expr = te_compile(funcao_str, vars, 1, &err);
    if (!expr) {
        printf("Erro na expressão, posição do erro: %d\n", err);
        return 0;
    }
    double val = te_eval(expr);
    te_free(expr);
    return val;
}

double derivada(double x) {
    double h = 1e-8;
    return (funcao(x + h) - funcao(x - h)) / (2*h);
}

double funcao_phi(double x) {
    return ((x*x*x)/9.0) + (1.0/3.0);
}

void bisseccao(double a, double b, double tol, int max_iter, FILE *fp) {
    fprintf(fp, "Método da Bissecção\n");
    fprintf(fp, "Iter\t   a\t\t   b\t\t   meio\t\t f(meio)\n");

    int k = 0;
    double meio;

    while (fabs(b - a) > tol && k < max_iter) {
        k++;
        double fa = funcao(a);
        double fb = funcao(b);
        meio = (a + b) / 2.0;
        double fmeio = funcao(meio);

        fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\n", k, a, b, meio, fmeio);

        if (fa * fmeio < 0)
		   	b = meio;
        else
		   	a = meio;
    }

    fprintf(fp, "k = %d\nRaiz = %lf\n\n", k, meio);
}

void iterativo_linear(double x0, double tol, int max_iter, FILE *fp) {
    fprintf(fp, "Método Iterativo Linear (Ponto Fixo)\n");
    fprintf(fp, "Iter\t   x_k\t\t   x_k+1\t\t   |x_k+1 - x_k|\n");

    int k = 0;
    double x1;

    while (k < max_iter) {
        x1 = funcao_phi(x0);
        fprintf(fp, "%d\t%lf\t%lf\t%lf\n", k + 1, x0, x1, fabs(x1 - x0));

        if (fabs(x1 - x0) < tol) {
            fprintf(fp, "k = %d\nRaiz aproximada: %lf\n\n", k + 1, x1);
            return;
        }

        x0 = x1;
        k++;
    }

    fprintf(fp, "Máximo de iterações atingido.\nÚltima aproximação: %lf\n\n", x1);
}

void newton_raphson(double x0, double tol, int max_iter, FILE *fp) {
    fprintf(fp, "Método de Newton-Raphson\n");
    fprintf(fp, "Iter\t     x\t\t    f(x)\t   |Δx|\n");

    int k = 0;
    double x1, fx0, dfx0;

    while (k < max_iter) {
        fx0 = funcao(x0);
        dfx0 = derivada(x0);

        if (fabs(fx0) < tol) {
            fprintf(fp, "Raiz aproximada encontrada: %lf\n", x0);
            return;
        }

        if (dfx0 == 0) {
            fprintf(fp, "Derivada igual a zero. Não é possível continuar.\n");
            return;
        }

        x1 = x0 - fx0 / dfx0;
        fprintf(fp, "%d\t%lf\t%lf\t%lf\n", k + 1, x0, fx0, fabs(x1 - x0));

        if (fabs(x1 - x0) < tol) {
            fprintf(fp, "k =  %d\nRaiz: %lf\n\n", k + 1, x1);
            return;
        }

        x0 = x1;
        k++;
    }

    fprintf(fp, "Máximo de iterações atingido.\nÚltima aproximação: %lf\n\n", x0);
}

void secante(double x0, double x1, double tol, int max_iter, FILE *fp) {
    fprintf(fp, "Método da Secante\n");
    fprintf(fp, "Iter\t    x0\t\t    x1\t\t   f(x1)\t    novo x\n");

    int iter = 0;
    double x_temp;

    while (iter < max_iter) {
        double fx0 = funcao(x0);
        double fx1 = funcao(x1);

        if (fabs(fx1) < tol) {
            fprintf(fp, "Raiz aproximada encontrada: %lf\n\n", x1);
            return;
        }

        if (fx1 == fx0) {
            fprintf(fp, "Erro: f(x1) = f(x0). Divisão por zero.\n");
            return;
        }

        x_temp = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\n", iter + 1, x0, x1, fx1, x_temp);

        if (fabs(x_temp - x1) < tol) {
            fprintf(fp, "k = %d\nRaiz = %lf\n\n", iter + 1, x_temp);
            return;
        }

        x0 = x1;
        x1 = x_temp;
        iter++;
    }

    fprintf(fp, "Máximo de iterações atingido.\nÚltima aproximação: %lf\n\n", x1);
}

void regula_falsi(double a, double b, double tol, int max_iter, FILE *fp) {
    fprintf(fp, "Método Regula Falsi (Falsa Posição)\n");
    fprintf(fp, "Iter\t    a\t\t    b\t\t    c\t\t   f(c)\n");

    int iter = 0;
    double c;

    double fa = funcao(a);
    double fb = funcao(b);

    if (fa * fb >= 0) {
        fprintf(fp, "A função deve mudar de sinal no intervalo dado.\n\n");
        return;
    }

    while (iter < max_iter) {
        c = b - fb * (b - a) / (fb - fa);
        double fc = funcao(c);

        fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\n", iter + 1, a, b, c, fc);

        if (fabs(fc) < 1e-8) {
            fprintf(fp, "Raiz = %lf\n\n", c);
            return;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }

        iter++;
    }

    fprintf(fp, "Máximo de iterações atingido.\nÚltima aproximação: %lf\n\n", c);
}

int main() {
    FILE *fp = fopen("resultados.txt", "w");
    if (!fp) {
        printf("Erro ao criar arquivo resultados.txt\n");
        return 1;
    }

    // Lê a função do arquivo funcao.txt
    FILE *finput = fopen("funcao.txt", "r");
    if (!finput) {
        printf("Erro ao abrir funcao.txt\n");
        return 1;
    }
    fgets(funcao_str, sizeof(funcao_str), finput);
    fclose(finput);

    printf("Função lida: %s\n", funcao_str);

    bisseccao(0.0, 1.0, 1e-5, 50, fp);
    iterativo_linear(0.5, 0.0005, 50, fp);
    newton_raphson(0.5, 1e-5, 50, fp);
    secante(0.0, 1.0, 0.0001, 50, fp);
    regula_falsi(0.0, 1.0, 0.0001, 50, fp);

    fclose(fp);
    printf("Resultados salvos em 'resultados.txt'\n");

    return 0;
}

