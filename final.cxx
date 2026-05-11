#include <stdio.h>
#include <math.h>

double f(double x)
{
	// (sin(x) + 1) / pow(x, 1.0/3.0);
	
	return 3.0 * x * (sin(x * x * x) + 1);
	
	//return 0.5* log(1.0 + x*x);
	
	//return atan(x*(x*x+1)/(sqrt(x*x*x*x+1)));
	
	//return 2.0 * (sin(x*x) + 1);
	
	//return log(pow(1 + pow(x,4), 1.0 / 3.0));
}

int sign(double x)
{
	if (x > 0) return 1;
	else if (x < 0) return -1;
	return 0;
}

double Si(double (*f)(double), double a, double b)
{
	return ((b - a) / 12.0) * (f(a) + 5 * f((a + b) / 2.0 - (b - a) / (2.0 * sqrt(5.0))) + 5 * f((a + b) / 2.0 + (b - a) / (2.0 * sqrt(5.0))) + f(b));
}

double Integrate2(double (*f)(double), double a, double b, double eps)
{
    double I = 0,  h = 0.1, s1, s2, delta, xi;
    
    while (a < b)
    {
        s1 = Si(f, a, a + h);
        s2 = Si(f, a, a + h / 2.0) + Si(f, a + h / 2.0, a + h);
        
        delta = (s2 - s1) / 63.0;

        if (fabs(delta) < eps)
        {
            if (a + h > b)
            {
            	h = b - a;
            	I += Si(f, a, a + h / 2.0) + Si(f, a + h / 2.0, a + h);
     	       a += h;
     	       return I;
            }
            
            I += s2;
            a += h;
        }
        
        xi = pow(fabs(delta)/eps, 1.0/7.0);
        if (xi > 10.0) xi = 10.0;
        if (xi < 0.1) xi = 0.1;
        
        h = 0.95 * h / xi;
    }
    
    return I;
}

double Integrate1(double (*f)(double), double *a, double A, double eps, double *Rg, double *Rs, int *kol, double *h)
{
	double I = 0, s1, s2, delta, xi, h_new = 0.1;
	
	while (I < A)
	{
		*h = h_new;

		s1 = Si(f, *a, *a + *h);
		s2 = Si(f, *a, *a + (*h) / 2.0) + Si(f, *a + (*h) / 2.0, *a + *h);

		delta = (s2 - s1) / 63.0;

		if (fabs(delta) < eps)
		{
			if (I + s2 > A) return I;
			
			I += s2;
			*Rg += fabs(delta);
			*Rs += delta;
			*a += *h;
			(*kol)++;
		}

		xi = pow(fabs(delta) / eps, 1.0 / 7.0);
		
		if (xi > 10.0) xi = 10.0;
		if (xi < 0.1) xi = 0.1;

		h_new = 0.95 * (*h) / xi;
	}

	return I;
}

double I2(double (*f)(double), double a, double h, double I, double A)
{
	return (I + Si(f, a, a + h  /  2.0) + Si(f, a + h  /  2.0, a + h)) - A;
}

int root_chords(double *x, double a, double h, double (*I2)(double (*)(double), double, double, double, double), double eps, double I, double A)
{
	double fa, fb, fc, c[2], b = a + h, a1 = a;
	int i = 0, ifa, ifb, ifc, k;

	c[0] = a; c[1] = b;
	fa = I - A; fb = I2(f, a, h, I, A);

	if (fabs(fa) < eps) {*x = a; return 1;}
	if (fabs(fb) < eps) {*x = a + h; return 1;}

	ifa = sign(fa); ifb = sign(fb);
	if (ifa * ifb != -1) return 0;

	for (k = 1; fabs(c[1] - c[0]) > eps; k++)
	{
		i = !i;

		c[i] = (a * fb - b * fa) / (fb - fa);
		fc = I2(f, a1, c[i] - a1, I, A);
		ifc = sign(fc);

		if (fabs(fc) < eps) {*x = c[i]; return k + 1;}

		if (ifc * ifa == 1) {a = c[i]; fa = fc;}

		else {b = c[i]; fb = fc;}
	}

	*x = c[i];

	return k;
}

int main()
{
	double R, Rg, Rs, I[3], eps[3] = {1e-7, 1e-9, 1e-11},  h, s1, s2, delta;
	double x_min, x_max, x, eps2 = 1e-12;

	double a = 0, A;
	int kolI = 0, kolX;
	
	scanf("%lf", &A);
	
	if(A < 0) {printf("alpha < 0"); return 1;}
	
	for(int i = 0; i < 3; i++)
	{
		Rs = 0; Rg = 0; kolI = 0; a = 0;
		
		I[i] = Integrate1(f, &a, A, eps[i], &Rg, &Rs, &kolI, &h);
	
		kolX = root_chords(&x, a, h, I2, eps2, I[i], A);
	
	
		s1 = Si(f, a, x);
		s2 = Si(f, a, (a + x) / 2.0) + Si(f, (a + x) / 2.0, x);
	
		delta = (s2 - s1) / 63.0;
	
		Rg += fabs(delta);		Rs += delta;		Rs = fabs(Rs);
	
		I[i] += s2;
	
		x_min = x - (I[i] - A) / f(x) - Rg / f(x);
		x_max = x - (I[i] - A) / f(x) + Rg / f(x);
		
		printf(" eps = %g:\n\n", eps[i]);
		printf(" Integral = %g\n x = %.6lf\n", I[i], x*x*x);
	
		printf(" x_min = %.6lf\n x_max = %.6lf\n", pow(x_min, 3), pow(x_max, 3));
	
		printf(" кол. итераций в методе хорд = %d\n кол. итераций в интегр. = %d\n",  kolX, kolI);
	
		printf(" Сред. погр. = %e\n Гаран. погр. = %e\n\n\n", Rs, Rg);
	
		//printf("%g		%g\n\n", pow(x, 3) - pow(x_min, 3), pow(x_max, 3) - pow(x, 3));
	}
	
	for(int i = 0; i < 3; i++)
	{
		I[i] = Integrate2(f, 0, x, eps[i]);
	}
	
	R = (I[0] - I[1]) / (I[1] - I[2]);
	
	printf(" коэф. сход. = %g\n", R);
	
	return 0;
}
