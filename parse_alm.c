#include <stdio.h> 
#include <math.h>
#include <time.h>

double get_double(FILE *alm_file){
    double var;
    fread( (void*)(&var), 8, 1, alm_file);
    return var;
}

int get_byte(FILE *alm_file){
    unsigned char var[1];
    fread( (void*)(&var), 1, 1, alm_file);
    int v2 = (int)var[0]; 
    return v2;
}

//СТРУКТУРА АЛЬМАНАХА ПОЛУЧЕННАЯ ИЗ ФАЙЛА
typedef struct gonec_alm {
        int probe;
        int coil;
        int day;
        int month;
        int year;
        double Te;
        double a0;  
        double e0;
        double i0;
        double omega0;
        double GDVU;
        double Tdr;
        double OMEGA_POINT; 
        double omega_point;        
        int days; 
} g_alm; 

g_alm get_alm_struct(FILE *alm_file){
    g_alm galm;
    //PROBE 
    galm.probe = get_byte(alm_file); 
    printf("probe:\t%d\n", galm.probe); 
    
    //COIL 
    unsigned char w[2]; 
    fread(w, 2, 1, alm_file);
    galm.coil = (w[0]<<8) + (int)w[1];
    unsigned char data[3];
    printf("coil:\t%d\n", galm.coil); 

    //DAY MONTH YEAR 
    fread(data, 3, 1, alm_file);
    galm.day = data[0];
    galm.month = data[1];
    galm.year = data[2];
    
    //TimeE OF CIRCLE OF ORIGEN 
    galm.Te = get_double(alm_file); 
    printf("Te: %f\n", galm.Te);           
    
    //
    galm.a0 = get_double(alm_file); 
    printf("a0: %f\n", galm.a0);           
    
    // 
    galm.e0 = get_double(alm_file); 
    printf("e0: %f\n", galm.e0);           
     
    // 
    galm.i0 = get_double(alm_file); 
    printf("div: %f\n", galm.a0);           
    
    galm.omega0 = get_double(alm_file); 
    printf("Perige: %f\n", galm.omega0);           
    
    galm.GDVU = get_double(alm_file); 
    printf("div: %f\n",galm.GDVU);           
    
    galm.Tdr = get_double(alm_file); 
    printf("Tdr: %f\n", galm.Tdr );           

    galm.OMEGA_POINT = get_double(alm_file); 
    printf("OMEGA_POINT: %f\n", galm.OMEGA_POINT);           
    
    galm.omega_point = get_double(alm_file); 
    printf("omega_point: %f\n", galm.omega_point);           
    
    galm.days = get_byte(alm_file);
    printf("days: %d\n", galm.days);    
    return galm;
}

time_t alm_to_time_t(int day, int month, int year, double TimeE){
    int hours = floor(TimeE / (60 * 60));
    int minutes = floor((TimeE - hours * 60 * 60) / 60);
    int seconds = TimeE - hours*60*60 - minutes*60;
    printf("HOURS: %d, MINUTES: %d, SECONDS: %d\n", hours, minutes, seconds);
    struct tm te;
    te.tm_sec = seconds;
    te.tm_min = minutes;
    te.tm_hour = hours;     
    te.tm_year = 2000 + year - 1900;
    te.tm_mon = month - 1;
    te.tm_mday = day;
    time_t e_time = mktime(&te);
    return e_time;
}
int main(){
    const int ALM_OFFSET = 11; 
    unsigned char b[1];
    int z = 0; 
    char * filename = "curalm.dat";
    
    FILE *alm_file;
    alm_file = fopen(filename, "rb");
    if(!alm_file){
        printf("%s", "file open error\n");
    }  
    else{
        printf("%s\n", "file successfully opened" );
    } 
    fseek(alm_file, ALM_OFFSET, SEEK_SET);
    g_alm galm; 
    galm = get_alm_struct(alm_file);
    //return 0; 
    
        
    //printf("TIME PASSED SECONDS: %d\ DAYS: %d \n", dif_time, days_passed); 
    /* 
    int  day = galm.day; 
    int  month = galm.month;
    int  year = galm.year;

    double Te = galm.Te;
    double a0 = galm.a0;
    double e0 = galm.e0;
    double i0 = galm.i0;
    double omega0 = galm.omega0; 
    double GDVU = galm.GDVU; 
    double Tdr = galm.Tdr;
    double OMEGA_POINT = galm.OMEGA_POINT;
    double omega_point = galm.omega_point; 
    int    days = galm.days; 
    */  
    
    int  day = 4; 
    int  month = 9;
    int  year = 16;
    double Te = 83756.689760;
    double a0 = 7875165.021004;
    double e0 = 0.001092;
    double i0 = 1.439208; 
    double omega0 = 1.162546; 
    double GDVU = 3.683670; 
    double Tdr = 6955.057288;
    double OMEGA_POINT = 7.304751/100000;
    double omega_point = 0.000478/100000; 
    
    
    //int  day = galm.day; 
    //int  month = galm.month;
    //int  year = galm.year;
    printf("Day %d Month %d Year %d\n", day, month, year); 
    
    time_t e_time = alm_to_time_t(day, month, year, Te); 
     
    //time_t c_time = time(NULL);
    
    

    time_t diff_time = 249 * Tdr;
    //int dif_time = c_time - e_time;
    //int days_passed =  floor(dif_time/(3600 * 24));

    int t = diff_time; 
    //int t = Tdr; 
    int t0 = 0; 
    double OMEGA = GDVU - OMEGA_POINT * (t - t0);
    double omega = omega0 + omega_point * (t - t0);
    double TRUE_ANOMALIA = -1 * omega;
    double E0 = 2 * atan( sqrt( (1 - e0) / (1 + e0) ) * tan(TRUE_ANOMALIA / 2) );
    double M0 = E0 - e0 * sin(E0);
    double n0 = (2 * M_PI) / Tdr;
    double M = M0 + n0 * (t - t0); 
    double Eprev = M;
    double Ecur = 0; 
    while(1){
        Ecur = M + e0 * sin(Eprev); 
        float delta_e = fabs(Ecur-Eprev); 
        printf("DELTA_E :\t%.12f\n", delta_e);
        if( delta_e >= 0.0000000001 ){
            printf("FINAL Ei:\t%f\n", Ecur);
            printf("FINAL Ei-1:\t%f\n", Eprev);
            Eprev = Ecur; 
            continue;  
        }  
        printf("FINAL Ei:\t%f\n", Ecur);
        printf("FINAL Ei-1:\t%f\n", Eprev);
        printf("DELTA :\t%.12f\n", delta_e);
        printf("%s", "break" );
        break;  
    }
    printf("FINAL E:\t%f\n", Ecur);
    
    double r = a0 * (1 - E0*cos(Ecur));
    double HUITA1, HUITA2;
    double Qx, Qy, Qz; 
    HUITA1 = a0 * (cos(Ecur) - e0);
    HUITA2 = a0 * sqrt(1 - pow(e0, 2) ) * sin(Ecur);
    double Px, Py, Pz; 
    Px = cos(omega) * cos(OMEGA) - sin(omega) * sin(OMEGA) * cos(i0);
    Py = cos(omega) * sin(OMEGA) + sin(omega) * cos(OMEGA) * cos(i0);
    Pz = sin(omega) * sin(i0);
    Qx = -1 * sin(omega) * cos(OMEGA) - cos(omega) * sin(OMEGA) * cos(i0);
    Qy = -1 * sin(omega) * sin(OMEGA) + cos(omega) * cos(OMEGA) * cos(i0);
    Qz = cos(omega) * sin(i0);
    
    double X = Px * HUITA1  +  Qx * HUITA2;
    double Y = Py * HUITA1  +  Qy * HUITA2;          
    double Z = Pz * HUITA1  +  Qz * HUITA2;           
    
    printf("==========================================\n");
    printf("X:\t%f\n", X); 
    printf("Y:\t%f\n", Y); 
    printf("Z:\t%f\n", Z); 
    printf("==========================================\n");
    double L = atan(Y/X);
    double alfa = 1.0/298.3; 
    double rr = sqrt( pow(X, 2) + pow(Y, 2) );  
    double B = atan( Z / ( pow((1 - alfa),2)*rr) );      
     printf("==========================================\n");
    printf("DOLGOTA: DEGREE\t%f\n", 360 * L/(2 * M_PI)); 
    printf("SHIROTA: DEGREE\t%f\n", 360 * B/(2 * M_PI)); 
    //printf("L: RAD   \t%f\n", L); 
    printf("==========================================\n");
    
    return 1;
}
