//Integral Masa+Energia+instante
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

const double pi = 3.14;
const double g = 9.81;

void Malla1D(double *x, const double L, const int N); //Cambiado
void Condicion_Inicial_Step(double *Matrix_Nk,double *x, const int N, const double L, const double h1, const double h2); //Condiciones inciales Nk=1

void SetBC_periodic(double *Xk, const int N);
void Calcular_RHS_U(double *rhs_u, double *Nk, const double g, const double Delta_X, const int N);
void Calcular_RHS_N(double *rhs_n, double *Uk, double *Nk, const double h, const double Delta_X, const int N);
void Avance_FowardEuler(double *Xkn1, double *Xk, double* rhs_u, double delta_t, const int N);
void ShapiroFilter(double *Xk, int N, double e);

void Imprimir(double *x, double *Uk, double *Nk, double t, int instant, int N);void Imprimir2(double *x, double *Uk, double *Nk, double t, int instant, int N);
void Imprimir3(double *x_h, double *x_u, double *Uk, double *Nk, double t, int instant, int N) {
	int i;
	char nombre[100];
	sprintf(nombre,"%d.txt",instant);
	
	FILE *f=fopen(nombre, "w");
	
	fprintf(f,"%f \n",t);	
	for(i = 0; i < N+2; ++i) //for(i=1;i<=N;++i)
		fprintf(f,"%f %f %f %f \n",x_h[i],x_u[i],Uk[i],Nk[i]);	
	
	fclose(f);
}


void ResetearMatrices(double *Ukn1, double *Nkn1, double *Uk, double *Nk,int N);

//Masa
double Integral(double *arr, int size,double base,double h) { // para conservacion masa
    int i;double area; double Nk1,Nk2;
  //  printf("%2.f ",base);
 
  	for (i = 1; i <= size; i++) {
      //  sum =sum + arr[i];
      	Nk1=arr[i]; Nk2=arr[i+1];
      //	printf("%f %f / %2.f \n",Nk1,Nk2,base);
    	area=area+(Nk1)*base+(Nk2-Nk1)*base/2;
    	if (Nk2>0){
    	//	printf(" error ");
		}
    }
    //printf("\n");
 //   area=(arr[i])*base+(arr[i+1]-arr[i])*base/2; 
  //  printf("%3.f \n",area);
    return area;
}


void ImprimirVC(double *vector, int length){
    FILE *file = fopen("IntregralEuler1.txt", "w"); // Open the file in write mode

    if (file == NULL) {
        printf("Error opening the file.\n");
        return;
    }
	int i;
    for ( i = 0; i < length; i++) {
        fprintf(file, "%lf\n", vector[i]); // Write each element of the vector to the file
    }

    fclose(file); // Close the file
}



//Energia
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

void ImprimirEn(double *vector, int length){
    FILE *file = fopen("EnergyEulerv2.txt", "w"); // Open the file in write mode

    if (file == NULL) {
        printf("Error opening the file.\n");
        return;
    }
	int i;
    for ( i = 0; i < length; i++) {
        fprintf(file, "%lf\n", vector[i]); // Write each element of the vector to the file
    }

    fclose(file); // Close the file
}

void ImprimirEn2(double *vector, int length,double *time){
    FILE *file = fopen("EnergyEulerv2.txt", "w"); // Open the file in write mode

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
double maxval(double *x, const int N){
	double m=-1e99; 
	int i;
	for( i=1; i<=N; ++i) m = fmax(m,x[i]); 
	return m; 
} //{ double m=-1e99; for(int i=1; i<=N; ++i) m = fmax(m,x[i]); return m; } 
double minval(double *x, const int N){
	double m=1e99;  
	int i;
	for( i=1; i<=N; ++i) m = fmin(m,x[i]); 
	return m;
} //{ double m=1e99;  for(int i=1; i<=N; ++i) m = fmin(m,x[i]); return m; }

int isapprox(double x1, double x2) { return fabs(x1 - x2) < 1e-6; }

int main(){	
	clock_t start, end;
    double cpu_time_used;

    start = clock(); // Start measuring time
	// Parametros
	const int    N       = 1000;    // numero de NODOS en X
	const double Long    = 1000;   // Longitud *metros
	const double Delta_X = Long/N; // Delta X
	printf("%f \n",Delta_X);
	
	const double Tiempo  = 200000.; // tiempo del experiemento // 100000
	const double T_save  = 0.5;  // AM: Cada cuanto queremos guardar 
	const double Delta_T = 0.1;  // AM: normalmente se fija el dt
	double dt=Delta_T;
	
	const double g       = 9.81; // gravedad (m*s^2) -gravedad
	const double e       = 0.1;   // Factor Disipacion
	
	const double h1      = 160.; // Anchura ola inicial  
	const double h2      = 1.;   // Altura ola inicial  
	const double h       = 10.;  // Altura nivel del mar -

	int inst = 0, isave = 1; // AM: instante
	double t = 0., t2 = T_save;
	
	//Matrices
	double *x, *Matrix_Uk, *Matrix_Ukn1, *Matrix_Nk, *Matrix_Nkn1, *rhs_u, *rhs_n;
	x           = (double*)calloc((N+2),sizeof(double)); 
	Matrix_Uk   = (double*)calloc((N+2),sizeof(double)); 
	Matrix_Nk   = (double*)calloc((N+2),sizeof(double));
	Matrix_Ukn1 = (double*)calloc((N+2),sizeof(double)); 
	Matrix_Nkn1 = (double*)calloc((N+2),sizeof(double)); 
	rhs_u       = (double*)calloc((N+2),sizeof(double));
	rhs_n       = (double*)calloc((N+2),sizeof(double)); 

	double *Masa; // conservacion masa
	
	Malla1D(x,Long,N);
	Condicion_Inicial_Step(Matrix_Nk,x,N,Long,h1,h2); 
/*	Masa = (double*)calloc((1),sizeof(double)); 
	Masa[inst]=Integral(Matrix_Nk,N,Delta_X);*/
	
	SetBC_periodic(Matrix_Nk,N);
	Imprimir(x,Matrix_Uk,Matrix_Nk,0.0,0,N);


	//Integral
	int cont2=0;
	
	double *integral;
	integral           = (double*)calloc((1),sizeof(double)); 
	
	//store times
	double *tiempos;
	tiempos           = (double*)calloc((1),sizeof(double)); 
	tiempos[0]=0;
	
	//Energy
	double *EnergyVect;
	EnergyVect           = (double*)calloc((1),sizeof(double)); 
	EnergyVect[0]=IntegralEnergy(Matrix_Uk,Matrix_Nk,h,N); //instante inicial
	
	while(1){
		inst++;	t = t + Delta_T;
		
		Calcular_RHS_U(rhs_u,Matrix_Nk,g,Delta_X,N);
		Avance_FowardEuler(Matrix_Ukn1,Matrix_Uk,rhs_u,Delta_T,N);
		SetBC_periodic(Matrix_Ukn1,N);
		
		Calcular_RHS_N(rhs_n,Matrix_Ukn1,Matrix_Nk,h,Delta_X,N);
		Avance_FowardEuler(Matrix_Nkn1,Matrix_Nk,rhs_n,Delta_T,N);
		ShapiroFilter(Matrix_Nkn1,N,e);
		SetBC_periodic(Matrix_Nkn1,N);
		
		int ins2=inst+1;
	//	Masa = (double*)realloc(Masa,ins2*sizeof(double));
	//	Masa[inst]=Integral(Matrix_Nkn1,N,Delta_X);
//		double integral= Integral(Matrix_Nkn1,N,Delta_X);
	//	printf("%f ",integral);
	
		integral=(double*)realloc(integral,(inst+1)*sizeof(double));
		integral[inst]=Integral(Matrix_Nkn1,N,Delta_X,h);
		
		EnergyVect=(double*)realloc(EnergyVect,(inst+1)*sizeof(double));
		EnergyVect[inst]=IntegralEnergy(Matrix_Ukn1,Matrix_Nkn1,h,N);

		//-----------Tiempo guardar
		tiempos=(double*)realloc(tiempos,(inst+1)*sizeof(double));
		tiempos[inst]=t;

	//	printf("%4i | dt = %.3e | t = %.3e | target = %.3e | umax = %.3e | umin = %.3e\n | Masa = %d ",inst,Delta_T,t,t2,maxval(Matrix_Ukn1,N),minval(Matrix_Ukn1,N),Masa[inst+1]);

		// AM: Guardamos un instante
		if ( isapprox(t,t2) ) {
	 	//	Imprimir2(x,Matrix_Ukn1,Matrix_Nkn1,t,isave,N);
	 //	Imprimir(x_h,x_u,Matrix_Ukn1,Matrix_Nkn1,t,isave,N);
			// Augmentar el tiempo de guardado
			t2 += T_save; isave++;
		} 
		
	
	//	printf("%f ",Masa[inst] );
		
		if (t >= Tiempo) break;
		
		
		
		ResetearMatrices(Matrix_Ukn1, Matrix_Nkn1, Matrix_Uk, Matrix_Nk, N);
	}
	int i;
	/*
	for (i = 0; i < inst; i++) {
    printf("%d ", Masa[i]);
	}
	*/	
//	ImprimirVC(integral,inst);
	ImprimirEn2(EnergyVect,inst,tiempos);
	
    end = clock(); // Stop measuring time
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Time elapsed during compilation: %f seconds\n", cpu_time_used);
	return 0;
}

void Malla1D(double *x, const double L, const int N){
	int i;
	double dx = L/N;
	x[0]   = -dx/2.;    // Por decir algo
	x[N+1] = L + dx/2.; // Por decir algo
	for(i = 1; i<= N; ++i){
		x[i] = x[i-1] + dx;
//		printf("%f ",x[i]);
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
		rhs_u[i] = g*(Nk[i+1] - Nk[i])/Delta_X; // Standard FD
	}
	
}

void Avance_FowardEuler(double *Xkn1, double *Xk, double* rhs_u, double delta_t, const int N) {
	int i;
	for (i = 1; i <= N; ++i){
		Xkn1[i] = Xk[i] - delta_t*rhs_u[i];
	}
}

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

		rhs_n[i] = (he - hw)/Delta_X;
	}
}

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

void ResetearMatrices(double *Ukn1, double *Nkn1, double *Uk, double *Nk,int N){
	int i;
	for(i = 0; i < N+2; ++i){
		Uk[i] = Ukn1[i];
		Nk[i] = Nkn1[i];
	}
}

void Imprimir(double *x, double *Uk, double *Nk, double t, int instant, int N) {
	int i;
	char nombre[100];
	sprintf(nombre,"%d.txt",instant);
	
	FILE *f=fopen(nombre, "w");
	
	fprintf(f,"%f \n",t);	
	for(i = 0; i < N+2; ++i) //for(i=1;i<=N;++i)
		fprintf(f,"%f %f %f \n",x[i],Uk[i],Nk[i]);	
	
	fclose(f);
}


void Imprimir2(double *x, double *Uk, double *Nk, double t, int instant, int N) {
    int i;
    char folder_name[] = "prueba2Euler-5T"; // Name of the folder to create
    char file_name[100];

    // Create folder if it doesn't exist
    #ifdef _WIN32
        // For Windows
        _mkdir(folder_name);
    #else
        // For Unix-like systems
        mkdir(folder_name, 0777);
    #endif

    sprintf(file_name, "%s/%d.txt", folder_name, instant); // Construct file path with folder name

    FILE *f = fopen(file_name, "w");
    if (f == NULL) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "%f\n", t);
    for (i = 0; i < N + 2; ++i)
        fprintf(f, "%f %f %f\n", x[i], Uk[i], Nk[i]);

    fclose(f);
}
