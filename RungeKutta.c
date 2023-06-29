//con integral
//con guardar integral
//con times y saved
//mas energy

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double pi = 3.14;
const double g = 9.81;

// Comment this define to switch to the old timestep
#define RKSTEPPING

#define RKMAXLEN 14
typedef struct _RKTable {
	int n_stages;
	double A[RKMAXLEN], B[RKMAXLEN], c[RKMAXLEN];
} RKTable;

RKTable rk4s5 = {
/*
Carpenter, M.H., Kennedy, C.A., 1994. 
Fourth-Order 2N-Storage Runge- Kutta Schemes (Technical Memorandum (TM) 
No. NASA-TM-109112). NASA.
*/
	.n_stages = 5,
	.A        = {0.0, -0.4178904745, -1.192151694643, -1.697784692471, -1.514183444257},
	.B        = {0.1496590219993, 0.3792103129999, 0.8229550293869, 0.6994504559488, 0.1530572479681},
	.c        = {0.0,0.1496590219993,0.3704009573644,0.6222557631345,0.9582821306748}
};

RKTable rk4s6 = {
/*
Allampalli, V., Hixon, R., Nallasamy, M., Sawyer, S.D., 2009. 
High-accuracy large-step explicit Runge–Kutta (HALE-RK) schemes for 
computational aeroacoustics. Journal of Computational Physics 228, 3837–3850. 
*/
	.n_stages = 6,
	.A        = {0.000000000000,-0.691750960670,-1.727127405211,-0.694890150986,-1.039942756197,-1.531977447611},
	.B        = {0.122000000000,0.477263056358,0.381941220320,0.447757195744,0.498614246822,0.186648570846},
	.c        = {0.000000000000,0.122000000000,0.269115878630,0.447717183551,0.749979795490,0.898555413085}
};
RKTable RK4s7 = {
	/*
	# Allampalli, V., Hixon, R., Nallasamy, M., Sawyer, S.D., 2009. 
			# High-accuracy large-step explicit Runge–Kutta (HALE-RK) schemes for 
			# computational aeroacoustics. Journal of Computational Physics 228, 3837–3850. 
			*/
			.n_stages = 7,
			.A = {.0,-0.647900745934,-2.704760863204,-0.460080550118,-0.500581787785,-1.906532255913,-1.450000000000},
			.B = {0.117322146869,0.503270262127,0.233663281658,0.283419634625,0.540367414023,0.371499414620,0.136670099385},
			.c = {.0,0.117322146869,0.294523230758,0.305658622131,0.582864148403,0.858664273599,0.868664273599}			
};
RKTable RK4s12 = {
	/*
	# Niegemann, J., Diehl, R., Busch, K., 2012. 
			# Efficient low-storage Runge–Kutta schemes with optimized stability regions. 
			# Journal of Computational Physics 231, 364–372.
	*/ 
	.n_stages = 12,
	.A = {0.0000000000000000,-0.0923311242368072,-0.9441056581158819,-4.3271273247576394,-2.1557771329026072,-0.9770727190189062,-0.7581835342571139,-1.7977525470825499,-2.6915667972700770,-4.6466798960268143,-0.1539613783825189,-0.5943293901830616},
	.B = {0.0650008435125904,0.0161459902249842,0.5758627178358159,0.1649758848361671,0.3934619494248182,0.0443509641602719,0.2074504268408778,0.6914247433015102,0.3766646883450449,0.0757190350155483,0.2027862031054088,0.2167029365631842},
	.c = {0.0000000000000000 ,0.0650008435125904 ,0.0796560563081853 ,0.1620416710085376 ,0.2248877362907778 ,0.2952293985641261 ,0.3318332506149405 ,0.4094724050198658 ,0.6356954475753369 ,0.6806551557645497 ,0.7143773712418350 ,0.9032588871651854}
};


RKTable RK4s13  = {
	/*
	# Niegemann, J., Diehl, R., Busch, K., 2012. 
			# Efficient low-storage Runge–Kutta schemes with optimized stability regions. 
			# Journal of Computational Physics 231, 364–372. 
			*/
			.n_stages = 13,
			.A = {0.0000000000000000,-0.6160178650170565,-0.4449487060774118,-1.0952033345276178,-1.2256030785959187,-0.2740182222332805,-0.0411952089052647,-0.1797084899153560,-1.1771530652064288,-0.4078831463120878,-0.8295636426191777,-4.7895970584252288,-0.6606671432964504},
			.B = {0.0271990297818803,0.1772488819905108,0.0378528418949694,0.6086431830142991,0.2154313974316100,0.2066152563885843,0.0415864076069797,0.0219891884310925,0.9893081222650993,0.0063199019859826,0.3749640721105318,1.6080235151003195,0.0961209123818189},
			.c = {0.0000000000000000,0.0271990297818803,0.0952594339119365,0.1266450286591127,0.1825883045699772,0.3737511439063931,0.5301279418422206,0.5704177433952291,0.5885784947099155,0.6160769826246714,0.6223252334314046,0.6897593128753419,0.9126827615920843}
		
	
};

RKTable RK4s14 ={
	/*
	# Niegemann, J., Diehl, R., Busch, K., 2012. 
			# Efficient low-storage Runge–Kutta schemes with optimized stability regions. 
			# Journal of Computational Physics 231, 364–372. 
			*/
			.n_stages = 14,
			.A = {0.0000000000000000,-0.7188012108672410,-0.7785331173421570,-0.0053282796654044,-0.8552979934029281,-3.9564138245774565,-1.5780575380587385,-2.0837094552574054,-0.7483334182761610,-0.7032861106563359,0.0013917096117681,-0.0932075369637460,-0.9514200470875948,-7.1151571693922548},
			.B = {0.0367762454319673,0.3136296607553959,0.1531848691869027,0.0030097086818182,0.3326293790646110,0.2440251405350864,0.3718879239592277,0.6204126221582444,0.1524043173028741,0.0760894927419266,0.0077604214040978,0.0024647284755382,0.0780348340049386,5.5059777270269628},
			.c = {0.0000000000000000 ,0.0367762454319673 ,0.1249685262725025 ,0.2446177702277698 ,0.2476149531070420 ,0.2969311120382472 ,0.3978149645802642 ,0.5270854589440328 ,0.6981269994175695 ,0.8190890835352128 ,0.8527059887098624 ,0.8604711817462826 ,0.8627060376969976 ,0.8734213127600976 }

};

RKTable lddrk2s5 = {
/*
Stanescu, D., Habashi, W.G., 1998. 
2N-Storage Low Dissipation and Dispersion Runge-Kutta Schemes for Computational Acoustics. 
Journal of Computational Physics 143, 674–681.
*/
	.n_stages = 5,
	.A        = {0.0,-0.6913065,-2.655155,-0.8147688,-0.6686587},
	.B        = {0.1,0.75,0.7,0.479313,0.310392},
	.c        = {0.0,0.1,0.3315201,0.4577796,0.8666528}
};

RKTable lddrk4s6 = {
/*
Stanescu, D., Habashi, W.G., 1998. 
2N-Storage Low Dissipation and Dispersion Runge-Kutta Schemes for Computational Acoustics. 
Journal of Computational Physics 143, 674–681.
*/
	.n_stages = 6,
	.A        = {0.0,-0.4919575,-0.8946264,-1.5526678,-3.4077973,-1.0742640},
	.B        = {0.1453095,0.4653797,0.4675397,0.7795279,0.3574327,0.15},
	.c        = {0.0,0.1453095,0.3817422,0.6367813,0.7560744,0.9271047}
};

void Malla1D(double *x, const double L, const int N); //Cambiado
void Condicion_Inicial_Step(double *Matrix_Nk,double *x, const int N, const double L, const double h1, const double h2); //Condiciones inciales Nk=1

void SetBC_periodic(double *Xk, const int N);
void Calcular_RHS_U(double *rhs_u, double *Nk, const double g, const double Delta_X, const int N);
void Calcular_RHS_N(double *rhs_n, double *Uk, double *Nk, const double h, const double Delta_X, const int N);
void Avance_FowardEuler(double *Xkn1, double *Xk, double* rhs_u, double delta_t, const int N);
void ShapiroFilter(double *Xk, int N, double e);

void Imprimir(double *x, double *Uk, double *Nk, double t, int instant, int N);
void Imprimir2(double *x, double *Uk, double *Nk, double t, int instant, int N);
void ResetearMatrices(double *Ukn1, double *Nkn1, double *Uk, double *Nk,int N);

double TimeStep(const double F, double *Uk, double *Nk, const double h, const double g, const double dx, const int N);
void RK(double *Ukn1, double *Nkn1, double *Uk, double *Nk, double *rhs_u, double *rhs_n, double g, double h0, double Delta_X, double Delta_T, const int N, RKTable *t);

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

//-------------
double Integral(double *arr, int size,double base,double h) { // para conservacion masa
    int i;double area; double Nk1,Nk2;
  //  printf("%2.f ",base);
 
  	for (i = 0; i <= size+1; i++) {
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
    FILE *file = fopen("IntregralRK1.txt", "w"); // Open the file in write mode

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
		double eta=Nk[i] ;double h=h0+Nk[i];
		
	//	Energy= Energy+0.5*h*(Uk[i]*Uk[i])+0.5*g*h*h+g*h*eta;
		Energy= Energy+0.5*h*(Uk[i]*Uk[i])+0.5*g*h*h+g*h0*h;	
	}
	return Energy;
	
}

void ImprimirEn(double *vector, int length){
    FILE *file = fopen("EnergyRKT10.txt", "w"); // Open the file in write mode

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
    FILE *file = fopen("Energylddrk4s6 -T3v2.txt", "w"); // Open the file in write mode

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

//Imprimir Time
void ImprimirTime(double *vector, int length){
    FILE *file = fopen("TiemposRK10.txt", "w"); // Open the file in write mode

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


//----------
int isapprox(double x1, double x2) { return fabs(x1 - x2) < 1e-6; }

int main(){	
	// Parametros
	const int    N       = 1000;     // numero de NODOS en X
	const double Long    = 1000;   // Longitud *metros
	const double Delta_X = Long/N; // Delta X
	
	const double Tiempo  = 100.; // tiempo del experiemento
	const double T_save  = 0.5;  // AM: Cada cuanto queremos guardar 
	double Delta_T;              // AM: normalmente se fija el dt ----- Cambiado
	
	const double g       = 9.81; // gravedad (m*s^2) -gravedad
	const double e       = 0.0;  // Factor Disipacion
	
	const double h1      = 160.; // Anchura ola inicial  
	const double h2      = 1.;   // Altura ola inicial  
	const double h       = 10.;  // Altura nivel del mar -

	int inst = 0, isave = 1; // AM: instante
	double t = 0., t2 = T_save;
	
	const double visc    = 0;
	const double factorT = 0.1; //Factor de integracion temporal delta T ---------------------------------------------------
	RKTable  rkt         = lddrk2s5;

	//Matrices
	double *x, *Matrix_Uk, *Matrix_Ukn1, *Matrix_Nk, *Matrix_Nkn1, *rhs_u, *rhs_n;
	x           = (double*)calloc((N+2),sizeof(double)); 
	Matrix_Uk   = (double*)calloc((N+2),sizeof(double)); 
	Matrix_Nk   = (double*)calloc((N+2),sizeof(double));
	Matrix_Ukn1 = (double*)calloc((N+2),sizeof(double)); 
	Matrix_Nkn1 = (double*)calloc((N+2),sizeof(double)); 
	rhs_u       = (double*)calloc((N+2),sizeof(double));
	rhs_n       = (double*)calloc((N+2),sizeof(double)); 
	
	//
	Malla1D(x,Long,N);
	Condicion_Inicial_Step(Matrix_Nk,x,N,Long,h1,h2); 
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
	EnergyVect[0]=IntegralEnergy(Matrix_Uk,Matrix_Nk,h,N);
	
	
	while(1){
		Delta_T =17*TimeStep(factorT,Matrix_Uk,Matrix_Nk,h,g,Delta_X,N);
		// Si el paso de tiempo es mayor al paso de tiempo de guardado
		// ajustamos el paso de tiempo al paso de tiempo de guardado
		if (t + Delta_T > t2) Delta_T = t2 - t;
		
		inst++;	t = t + Delta_T;cont2++;

		#ifdef RKSTEPPING
		RK(Matrix_Ukn1,Matrix_Nkn1,Matrix_Uk,Matrix_Nk,rhs_u,rhs_n,g,h,Delta_X,Delta_T,N,&rkt);
		#else
		Calcular_RHS_U(rhs_u,Matrix_Nk,g,Delta_X,N);
		Avance_FowardEuler(Matrix_Ukn1,Matrix_Uk,rhs_u,Delta_T,N);
		SetBC_periodic(Matrix_Ukn1,N);
		
		Calcular_RHS_N(rhs_n,Matrix_Ukn1,Matrix_Nk,h,Delta_X,N);
		Avance_FowardEuler(Matrix_Nkn1,Matrix_Nk,rhs_n,Delta_T,N);
		SetBC_periodic(Matrix_Nkn1,N);
		# endif
		ShapiroFilter(Matrix_Nkn1,N,e);
		
		printf("%4i | dt = %.3e | t = %.3e | target = %.3e | umax = %.3e | umin = %.3e\n",inst,Delta_T,t,t2,maxval(Matrix_Ukn1,N),minval(Matrix_Ukn1,N));
/*
		integral=(double*)realloc(integral,(inst+1)*sizeof(double));
		integral[inst]=Integral(Matrix_Nkn1,N,Delta_X,h);
*/		
		EnergyVect=(double*)realloc(EnergyVect,(inst+1)*sizeof(double));
		EnergyVect[inst]=IntegralEnergy(Matrix_Ukn1,Matrix_Nkn1,h,N);
	//	printf("%.2f ",integral);
	
	//-----------Tiempo guardar
		tiempos=(double*)realloc(tiempos,(inst+1)*sizeof(double));
		tiempos[inst]=t;
	//	printf("%f ",t);
		
		// AM: Guardamos un instante
		if ( isapprox(t,t2) ) {
	 //		Imprimir2(x,Matrix_Ukn1,Matrix_Nkn1,t,isave,N);
			// Augmentar el tiempo de guardado
			t2 += T_save; isave++;
		} 
		//cada x instante
		/*
		if(cont2==3){
			double integral= Integral(Matrix_Nkn1,N,Delta_X, h);
		printf("%.2f \n",integral);
			cont2=0;
		}
		*/
		//al instante
		/*
		double integral= Integral(Matrix_Nkn1,N,Delta_X, h);
		printf("%.2f \n",integral);
		*/
		if (t >= Tiempo) break;
		
		ResetearMatrices(Matrix_Ukn1, Matrix_Nkn1, Matrix_Uk, Matrix_Nk, N);
	//	Imprimir2(x,Matrix_Ukn1, Matrix_Nkn1,t,inst,N);//
	}
//	ImprimirVC(integral,inst);
	ImprimirEn2(EnergyVect,inst,tiempos);
//	ImprimirTime(tiempos,inst);
	printf("%d ",inst);
	return 0;
}

void Malla1D(double *x, const double L, const int N){
	int i;
	double dx = L/N;
	x[0]   = -dx/2.;    // Por decir algo
	x[N+1] = L + dx/2.; // Por decir algo
	for(i = 1; i<= N; ++i)
		x[i] = x[i-1] + dx;
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

void Avance_FowardEuler(double *Xkn1, double *Xk, double* rhs_u, double delta_t, const int N) {
	int i;
	for (i = 1; i <= N; ++i){
		Xkn1[i] = Xk[i] + delta_t*rhs_u[i];
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

		rhs_n[i] = -(he - hw)/Delta_X;
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

//
void RK(double *Ukn1, double *Nkn1, double *Uk, double *Nk, double *rhs_u, double *rhs_n, double g, double h, double Delta_X, double dt, const int N, RKTable *t) {
	// Using Ukn1,Uk as Ku1,Ku2
	// Using Nkn1,Uk as Kn1,Kn2

	int i, s;

	// Copy Uk,Nk to Ukn1,Nkn1 
	for(i=0;i<N+2;++i) {
		Ukn1[i] = Uk[i]; // Ku1 = Uk
		Nkn1[i] = Nk[i]; // Kn1 = Nk
		// Set Uk and Nk to zero
		Uk[i]   = 0.;    // Ku2 = 0.
		Nk[i]   = 0.;    // Kn2 = 0.
	}

	// Loop number of stages
	for(s=0;s<t->n_stages;++s) {
		// 1. Compute the RHS of U and N equations
		Calcular_RHS_U(rhs_u,Nkn1,g,Delta_X,N);
		Calcular_RHS_N(rhs_n,Ukn1,Nkn1,h,Delta_X,N);
		// 2. Advance U and N to the next stage
		for(i=1;i<=N;++i) {
			Uk[i]    = t->A[s]*Uk[i] + dt*rhs_u[i]; // Ku2 = A*Ku2 + dt*rhs_u(Ku1)
			Ukn1[i] += t->B[s]*Uk[i];               // Ku1 = Ku1 + B*Ku2
			Nk[i]    = t->A[s]*Nk[i] + dt*rhs_n[i]; // Kn2 = A*Kn2 + dt*rhs_n(Kn1)
			Nkn1[i] += t->B[s]*Nk[i];               // Kn1 = Kn1 + B*Kn2
		}
		// 3. Enforce correct boundary conditions
		SetBC_periodic(Ukn1,N); // Now Ukn1 contains Ku1 for the next stage
		SetBC_periodic(Nkn1,N); // Now Nkn1 contains Kn1 for the next stage
	}
	
	
	// At the end of the function Uk and Nk are overwritten per stage
	// while Ukn1 and Nkn1 contain the next stage values
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

void Imprimir2(double *x, double *Uk, double *Nk, double t, int instant, int N) {
    int i;
    char folder_name[] = "lddrk4s6-10T"; // Name of the folder to create
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




