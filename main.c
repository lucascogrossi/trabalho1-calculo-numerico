#include <stdio.h>
#include <math.h>

double funcao(double x) {
	 return x*x*x - 9.0*x + 3.0;
}

double derivada(double x) {
    return 3.0*x*x - 9.0;
}

double funcao_phi(double x) {
	return ((x*x*x)/9.0) + (1.0/3.0);
}

void regula_falsi(double a, double b, double tol, int max_iter) {
    printf("Método Regula Falsi (Falsa Posição)\n");
    printf("Iter\t    a\t\t    b\t\t    c\t\t   f(c)\n");

    double fa = funcao(a);
    double fb = funcao(b);

    if (fa * fb >= 0) {
        printf("A função deve mudar de sinal no intervalo dado.\n\n");
        return;
    }

    int iter = 0;
    double c, fc;

    while (iter < max_iter) {
        c = b - fb * (b - a) / (fb - fa);
        fc = funcao(c);

        printf("%d\t%lf\t%lf\t%lf\t%lf\n", iter + 1, a, b, c, fc);

        if (fabs(fc) < tol) {
            printf("Raiz = %lf\n\n", c);
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

    printf("Máximo de iterações atingido.\n");
    printf("Última aproximação: %lf\n\n", c);
}


void secante(double x0, double x1, double tol, int max_iter) {
		printf("Método da Secante\n");
		printf("Iter\t    x0\t\t    x1\t\t   f(x1)\t    novo x\n");

		int iter = 0;
		double x_temp;

		while (iter < max_iter) {
			double fx0 = funcao(x0);
			double fx1 = funcao(x1);

			if (fabs(fx1) < tol) {
				printf("Raiz aproximada encontrada: %lf\n\n", x1);
				return;
			}

			if (fx1 == fx0) {
				printf("Erro: f(x1) = f(x0). Divisão por zero.\n");
				return;
			}

        x_temp = x1 - fx1 * (x1 - x0) / (fx1 - fx0);

        printf("%d\t%lf\t%lf\t%lf\t%lf\n", iter + 1, x0, x1, fx1, x_temp);

        if (fabs(x_temp - x1) < tol) {
            printf("k = %d\n", iter + 1);
            printf("Raiz = %lf\n\n", x_temp);
            return;
        }

        x0 = x1;
        x1 = x_temp;
        iter++;
    }

    printf("Máximo de iterações atingido.\n");
    printf("Última aproximação: %lf\n\n", x1);
}

void newton_raphson(double x0, double tol, int max_iter) {
    printf("Método de Newton-Raphson\n");
    printf("Iter\t     x\t\t    f(x)\t   |Δx|\n");

    int k = 0;
    double x1, fx0, dfx0;

    while (k < max_iter) {
        fx0 = funcao(x0);
        dfx0 = derivada(x0);

        if (fabs(fx0) < tol) {
            printf("Raiz aproximada encontrada: %lf\n", x0);
            return;
        }
        if (dfx0 == 0) {
            printf("Derivada igual a zero. Não é possível continuar.\n");
            return;
        }

        x1 = x0 - fx0 / dfx0;

        printf("%d\t%lf\t%lf\t%lf\n", k + 1, x0, fx0, fabs(x1 - x0));

        if (fabs(x1 - x0) < tol) {
            printf("k =  %d\n", k + 1);
            printf("Raiz: %lf\n\n", x1);
            return;
        }

        x0 = x1;
        k++;
    }

    printf("Máximo de iterações atingido.\n");
    printf("Última aproximação: %lf\n\n", x0);
}

void iterativo_linear(double x0, double tol, int max_iter) {
    printf("Método Iterativo Linear (Ponto Fixo)\n");
    printf("Iter\t   x_k\t\t   x_k+1\t\t   |x_k+1 - x_k|\n");

    int k = 0;
    double x1;

    while (k < max_iter) {
        x1 = funcao_phi(x0);
        printf("%d\t%lf\t%lf\t%lf\n", k + 1, x0, x1, fabs(x1 - x0));

        if (fabs(x1 - x0) < tol) {
            printf("k = %d\n", k + 1);
            printf("Raiz aproximada: %lf\n\n", x1);
            return;
        }

        x0 = x1;
        k++;
    }

    printf("Máximo de iterações atingido.\n");
    printf("Última aproximação: %lf\n\n", x1);
}



void bisseccao(double a, double b, double tol, int max_iter) {
	printf("Método da Bissecção\n");
	printf("Iter\t   a\t\t   b\t\t   meio\t\t f(meio)\n");
    int k = 0;
    if (fabs(b - a) < tol) {
        printf("Raiz = %lf\n", a);
    } else {
        double meio;
        while (fabs(b - a) > tol && k < max_iter) {
            k++;
            double finicio = funcao(a);
            meio = (a + b) / 2.0;
            double fmeio = funcao(meio);
			
			printf("%d\t%lf\t%lf\t%lf\t%lf\n", k, a, b, meio, fmeio);
            if (finicio * fmeio < 0) {
                b = meio;
            } else {
                a = meio;
            }
        }
        printf("k = %d\n", k);
        printf("Raiz = %lf\n\n", meio);
    }
}

int main() {
	// parametros bissecao
    double a_bi = 0.0;
    double b_bi = 1.0;
	double delta_bi = 1e-5;
    int max_iter_bi = 50;
	
    bisseccao(a_bi, b_bi, delta_bi, max_iter_bi);

	// parametros iterativo linear
	double x0_il = 0.5;
	double delta_il = 0.0005;
		int n_il = 50;

	iterativo_linear(x0_il, delta_il, n_il);

	// parametros metodo de newton
	double x0_newton = 0.5;
    double tol_newton = 1e-5;
    int max_iter_newton = 50;

    newton_raphson(x0_newton, tol_newton, max_iter_newton);	
	
	// parametros metodo secante
	double x0_sec = 0.0;
    double x1_sec = 1.0;
    double tol_sec = 0.0001;
    int max_iter_sec = 50;

	secante(x0_sec, x1_sec, tol_sec, max_iter_sec);
	
	// parametros metodo regula falsi
	double a_rf = 0.0;
	double b_rf = 1.0;
	double tol_rf = 0.0001;
	int max_iter_rf = 50;

	regula_falsi(a_rf, b_rf, tol_rf, max_iter_rf);

    return 0;
}

