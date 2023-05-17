#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>     /* setenv, getenv */
#include <string.h>

void c_setgetenv(char *emin, char *emax)
{
    //Set up all variable names:
    char* p;
    bool verbose = 1;
    char emin_string[30],emax_string[30];
  
    setenv("MU_ZONES","1",1);
    p = getenv("MU_ZONES");
    
    if (verbose == 1){
        if (p!=NULL){
            printf (" MU_ZONES is: %s\n",p);
        } else {
            printf (" Error setting MU_ZONES");
        }
    }
    
    setenv("ION_ZONES","20",1);
    p = getenv("ION_ZONES");
    
    if (verbose == 1){
        if (p!=NULL){
            printf (" ION_ZONES is: %s\n",p);
        } else {
            printf (" Error setting ION_ZONES");
        }
    }
        
    setenv("A_DENSITY","0",1);
    p = getenv("A_DENSITY");
    
    if (verbose == 1){
        if (p!=NULL){
            printf (" A_DENSITY is: %s\n",p);
        } else {
            printf (" Error setting A_DENSITY");
        }
    }
    
    setenv("REV_VERB","2",1);
    p = getenv("REV_VERB");
    
    if (verbose == 1){
        if (p!=NULL){
            printf (" REV_VERB is: %s\n",p);
        } else {
            printf (" Error setting REV_VERB");
        }
    }
        
    setenv("REF_VAR","1",1);
    p = getenv("REF_VAR");
    
    if (verbose == 1){
        if (p!=NULL){
            printf (" REF_VAR is: %s\n",p);
        } else {
            printf (" Error setting REF_VAR");
        }
    }
    
    setenv("ION_VAR","1",1);
    p = getenv("ION_VAR");
    
    if (verbose == 1){    
        if (p!=NULL){
            printf (" ION_VAR is: %s\n",p);
        } else {
            printf (" Error setting ION_VAR");
        }
    }
    
    setenv("EMIN_REF",emin,1);
    p = getenv("EMIN_REF");

    if (verbose == 1){    
        if (p!=NULL){
            printf (" EMIN_REF is: %s\n",p);
        } else {
            printf (" Error setting EMIN_REF");
        }
    }
    
    setenv("EMAX_REF",emax,1);
    p = getenv("EMAX_REF");
    
    if (verbose == 1){    
        if (p!=NULL){
            printf (" EMAX_REF is: %s\n",p);
        } else {
            printf (" Error setting EMAX_REF");
        }
    }
}
