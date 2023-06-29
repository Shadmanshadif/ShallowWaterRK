#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double g=9.81;
// Function prototypes
void Malla1D(double *x_h, double *x_u, const double L, const int N);
void Condicion_Inicial_Step(double *Matrix_Nk,double *x, const int N, const double L, const double h1, const double h2); //Condiciones inciales Nk=1
void SetBC_periodic(double *Xk, const int N);

void Calcular_RHS_U(double *rhs_u, double *Nk, const double g, const double Delta_X, const int N);
void Calcular_RHS_N(double *rhs_n, double *Uk, double *Nk, const double h, const double Delta_X, const int N);

void Avance_FowardEuler(double *Xkn1, double *Xk, double* rhs_u, double delta_t, const int N);
void AdamsBashforth2(double *Xkn1, double *Xk, double* rhs, double *rhs1, double delta_t, const int first, const int N);
void AdamsBashforth3(double *Xkn1, double *Xk, double* rhs, double *rhs1, double *rhs2, double delta_t, const int first, const int second, const int N);

double TimeStep(const double F, double *Uk, double *Nk, const double h, const double g, const double dx, const int N);

void ShapiroFilter(double *Xk, int N, double e);
void Imprimir(double *x_h, double *x_u, double *Uk, double *Nk, double t, int instant, int N);
void ResetearMatrices(double *Ukn1, double *Nkn1, double *Uk, double *Nk,int N);

double maxval(double *x, const int N){
	double m=-1e99; 
	int i;
	for( i=1; i<=N; ++i) m = fmax(m,x[i]); 
	return m; 
} 
double minval(double *x, const int N){
	double m=1e99;  
	int i;
	for( i=1; i<=N; ++i) m = fmin(m,x[i]); 
	return m;
} 

int isapprox(double x1, double x2) { return fabs(x1 - x2) < 1e-6; }

//Energy

double IntegralEnergy(double *Uk,double *Nk, double h0,int N){
	int i; double Energy=0;  
	
	for (i=0;i<=N+1;i++){
		double Nm=(Nk[i+1]+Nk[i])/2;
		double eta=Nm ;double h=h0+Nm;

		
		//Energy= Energy+0.5*h*(Uk[i]*Uk[i])+0.5*g*h0*h+g*h*eta;
		Energy= Energy+0.5*h*(Uk[i]*Uk[i])+0.5*g*h*h+g*h0*h;	
	}
	return Energy;
	
}

void ImprimirEn2(double *vector, int length,double *time){
    FILE *file = fopen("EnergyAB-T2v2.txt", "w"); // Open the file in write mode

    if (file == NULL) {
        printf("Error opening the file.\n");
        return;
    }
	int i;
    for ( i = 0; i < length; i++) {
        fprintf(file, "%lf  %lf\n",time[i], vector[i]); // Write each element of the vector to the file
    }

    fclose(file); // Close the file
}


//


int main(){	
	// Parametros
	const int    N       = 1000;     // numero de NODOS en X
	const double Long    = 1000;   // Longitud *metros
	const double Delta_X = Long/N; // Delta X
	
	const double Tiempo  = 200000.;   // tiempo del experiemento
	const double T_save  = 0.5;    // Cada cuanto queremos guardar 
	
	const double g       = 9.81;   // gravedad (m*s^2) -gravedad
	const double e       = 0.0;    // Factor Disipacion
	
	const double h1      = 160.;   // Anchura ola inicial  
	const double h2      = 1.;     // Altura ola inicial  
	const double h       = 10.;    // Altura nivel del mar -
	
	const double visc    = 0;
	const double factorT = 0.1;    // Factor de integracion temporal delta T

	
	// Variables auxiliares
	int inst = 0, isave = 1;
	double t = 0., t2 = T_save, Delta_T;

	// Matrices
	double *x_h, *x_u, *Matrix_Uk, *Matrix_Ukn1, *Matrix_Nk, *Matrix_Nkn1;
	double *rhs_u, *rhs_u1, *rhs_u2, *rhs_n, *rhs_n1, *rhs_n2;
	x_h         = (double*)calloc((N+2),sizeof(double)); 
	x_u         = (double*)calloc((N+2),sizeof(double)); 
	Matrix_Uk   = (double*)calloc((N+2),sizeof(double)); 
	Matrix_Nk   = (double*)calloc((N+2),sizeof(double));
	Matrix_Ukn1 = (double*)calloc((N+2),sizeof(double)); 
	Matrix_Nkn1 = (double*)calloc((N+2),sizeof(double)); 
	rhs_u       = (double*)calloc((N+2),sizeof(double));
	rhs_n       = (double*)calloc((N+2),sizeof(double)); 

	rhs_u1      = (double*)calloc((N+2),sizeof(double));
	rhs_u2      = (double*)calloc((N+2),sizeof(double));
	rhs_n1      = (double*)calloc((N+2),sizeof(double)); 
	rhs_n2      = (double*)calloc((N+2),sizeof(double));

	
	//
	Malla1D(x_h,x_u,Long,N);
	Condicion_Inicial_Step(Matrix_Nk,x_h,N,Long,h1,h2); 
	SetBC_periodic(Matrix_Nk,N);
	Imprimir(x_h,x_u,Matrix_Uk,Matrix_Nk,0.0,0,N);
	
	
	//Energy
	double *EnergyVect;
	EnergyVect           = (double*)calloc((1),sizeof(double));
	EnergyVect[0]=IntegralEnergy(Matrix_Uk,Matrix_Nk,h,N);  //
	//store times
	double *tiempos;
	tiempos           = (double*)calloc((1),sizeof(double)); 
	tiempos[0]=0;
	
	while(1){
		Delta_T = 2*TimeStep(factorT,Matrix_Uk,Matrix_Nk,h,g,Delta_X,N);
		// Si el paso de tiempo es mayor al paso de tiempo de guardado
		// ajustamos el paso de tiempo al paso de tiempo de guardado
		if (t + Delta_T > t2) Delta_T = t2 - t;
		inst++;	t = t + Delta_T;

		
		Calcular_RHS_U(rhs_u,Matrix_Nk,g,Delta_X,N);
		AdamsBashforth3(Matrix_Ukn1,Matrix_Uk,rhs_u,rhs_u1,rhs_u2,Delta_T,inst==1,inst==2,N);
		SetBC_periodic(Matrix_Ukn1,N);
		
		Calcular_RHS_N(rhs_n,Matrix_Ukn1,Matrix_Nk,h,Delta_X,N);
		AdamsBashforth3(Matrix_Nkn1,Matrix_Nk,rhs_n,rhs_n1,rhs_n2,Delta_T,inst==1,inst==2,N);
		SetBC_periodic(Matrix_Nkn1,N);
		
		ShapiroFilter(Matrix_Nkn1,N,e);
		
		printf("%4i | dt = %.3e | t = %.3e | target = %.3e | umax = %.3e | umin = %.3e\n",inst,Delta_T,t,t2,maxval(Matrix_Ukn1,N),minval(Matrix_Ukn1,N));
		
	

		EnergyVect=(double*)realloc(EnergyVect,(inst+1)*sizeof(double));
		EnergyVect[inst]=IntegralEnergy(Matrix_Ukn1,Matrix_Nkn1,h,N);
		
		//-----------Tiempo guardar
		tiempos=(double*)realloc(tiempos,(inst+1)*sizeof(double));
		tiempos[inst]=t;
		
		// AM: Guardamos un instante
		if ( isapprox(t,t2) ) {
	// 		Imprimir(x_h,x_u,Matrix_Ukn1,Matrix_Nkn1,t,isave,N);
			// Augmentar el tiempo de guardado
			t2 += T_save; isave++;
		} 
		if (t >= Tiempo) break;
		
		// Avanzar siguiente instante
		
		memcpy(rhs_n2,rhs_n1,(N+2)*sizeof(double));
		memcpy(rhs_n2,rhs_n1,(N+2)*sizeof(double));
		memcpy(rhs_u1,rhs_u,(N+2)*sizeof(double));
		memcpy(rhs_n1,rhs_n,(N+2)*sizeof(double));
	
		ResetearMatrices(Matrix_Ukn1, Matrix_Nkn1, Matrix_Uk, Matrix_Nk, N);
	}
	
		ImprimirEn2(EnergyVect,inst,tiempos);

	free(x_h); free(x_u); 
	free(Matrix_Uk); free(Matrix_Ukn1); 
	free(Matrix_Nk); free(Matrix_Nkn1);
	free(rhs_u); free(rhs_n); 
	#ifdef ABSTEPPING
	free(rhs_u1); free(rhs_u2); 
	free(rhs_n1); free(rhs_n2);
	#endif

	return 0;
}

void Malla1D(double *x_h, double *x_u, const double L, const int N){
	int i;
	double dx = L/N;
	x_h[0]    = -dx/2.;    // Por decir algo
	x_h[N+1]  = L + dx/2.; // Por decir algo
	x_u[0]    = 0.;        // Por decir algo
	x_u[N+1]  = L + dx;    // Por decir algo
	for(i = 1; i<= N; ++i) {
		x_h[i] = x_h[i-1] + dx;
		x_u[i] = x_u[i-1] + dx;
	}
}

void Condicion_Inicial_Step(double *Matrix_Nk, double *x, const int N, const double L, const double h1, const double h2){
	int i;
	for (i = 1; i <= N; ++i) { // Recuerda que 0 y N+1 son para BCs
		if (x[i] >= 0.5*(L-h1) && x[i] < 0.5*(L+h1)) {
			Matrix_Nk[i] += h2;
		}
	}
}
	
void SetBC_periodic(double *Xk, const int N){
	Xk[0]   = Xk[N];
	Xk[N+1] = Xk[1];
	
}

void Calcular_RHS_U(double *rhs_u, double *Nk, const double g, const double Delta_X, const int N){
	int i,k;
	// Loop en la malla
	for (i = 1; i <= N; ++i){ 
		rhs_u[i] = -g*(Nk[i+1] - Nk[i])/Delta_X; // Standard FD
	}
	
}

//
void Calcular_RHS_N(double *rhs_n, double *Uk, double *Nk, const double h0, const double Delta_X, const int N){
	int i;
	double hep, hen, he, hwp, hwn, hw;
	// Loop en la malla
	for (i = 1; i <= N; ++i) { 
		// depth grid points of a control volume, previous and next time step and
		// east and west sides of the control volume cell, pg. 72 eq. 4.19
		hep = 0.5*(Uk[i]   + fabs(Uk[i]))*(Nk[i]   + h0);   // previous 
		hen = 0.5*(Uk[i]   - fabs(Uk[i]))*(Nk[i+1] + h0);   // next time step
		he  = hep + hen; 								   // east
		hwp = 0.5*(Uk[i-1] + fabs(Uk[i-1]))*(Nk[i-1] + h0); // previous
		hwn = 0.5*(Uk[i-1] - fabs(Uk[i-1]))*(Nk[i]   + h0); // next 
		hw  = hwp + hwn;                                    // west

		rhs_n[i] = -(he - hw)/Delta_X;
	}
}

//
void Avance_FowardEuler(double *Xkn1, double *Xk, double* rhs, double delta_t, const int N) {
	int i;
	for (i = 1; i <= N; ++i){
		Xkn1[i] = Xk[i] + delta_t*rhs[i];
	}
}

void AdamsBashforth2(double *Xkn1, double *Xk, double* rhs, double *rhs1, double delta_t, const int first, const int N) {
	int i;
	if (first)
		// rhs1 is not available
		for (i = 1; i <= N; ++i) {
			Xkn1[i] = Xk[i] + delta_t*(1.5*rhs[i] - 0.5*rhs[i]);
		}
	else
		// rhs1 is available
		for (i = 1; i <= N; ++i) {
			Xkn1[i] = Xk[i] + delta_t*(1.5*rhs[i] - 0.5*rhs1[i]);
		}
}

void AdamsBashforth3(double *Xkn1, double *Xk, double* rhs, double *rhs1, double *rhs2, double delta_t, const int first, const int second, const int N) {
	int i;
	if (second)
		// rhs1 and rhs2 are not available, try with AB2
		AdamsBashforth2(Xkn1,Xk,rhs,rhs1,delta_t,first,N);
	else
		// rhs1 and rhs2 are available
		for (i = 1; i <= N; ++i) {
			Xkn1[i] = Xk[i] + delta_t*(23./12.*rhs[i] - 4./3.*rhs1[i] + 5./12.*rhs2[i]);
		}
}



double TimeStep(const double F, double *Uk, double *Nk, const double h0, const double g, const double dx, const int N) {
	double dtmin = 1e9;
	// Calcular maximos
	double Umax = maxval(Uk,N); // Alerta debe ser N por como esta programado maxval, sino hay overflow!
	double Hmax = maxval(Nk,N) + h0;
	// Calcular el paso mas restrictivo
	dtmin = fmin(dtmin, (isapprox(Umax,0.)) ? 1e9 : dx/Umax);
	dtmin = fmin(dtmin, dx/sqrt(g*Hmax));
	// Devolver el paso de tiempo multiplicado por el factor de seguridad
	return F*dtmin;
}

//
void ShapiroFilter(double *Xk, int N, double e) {
	if (isapprox(e,0.)) return; // No calcular si e == 0.
	int i;
	double buf[3];
	// Cargar datos en el buffer
	buf[0] = Xk[0];
	buf[1] = Xk[1];
	buf[2] = Xk[2];
	// Loop en la malla
	for (i = 1; i <= N; ++i) { 
		// Calcular Shapiro con el buffer
		Xk[i] = (1 - e)*buf[1] + 0.5*e*(buf[0] + buf[2]);
		// Update del buffer
		buf[0] = buf[1];
		buf[1] = buf[2];
		buf[2] = Xk[i+2];
	}
}

//
void ResetearMatrices(double *Ukn1, double *Nkn1, double *Uk, double *Nk,int N){
	int i;
	for(i = 0; i < N+2; ++i){
		Uk[i] = Ukn1[i];
		Nk[i] = Nkn1[i];
	}
}

void Imprimir(double *x_h, double *x_u, double *Uk, double *Nk, double t, int instant, int N) {
	int i;
	char nombre[100];
	sprintf(nombre,"%d.txt",instant);
	
	FILE *f=fopen(nombre, "w");
	
	fprintf(f,"%f \n",t);	
	for(i = 0; i < N+2; ++i) //for(i=1;i<=N;++i)
		fprintf(f,"%f %f %f %f \n",x_h[i],x_u[i],Uk[i],Nk[i]);	
	
	fclose(f);
}
