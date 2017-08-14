#include <iostream>


#include "~/mosek/7/tools/platform/linux64x86/h/mosek.h"
#include "~/mosek/7/tools/platform/linux64x86/bin/libmosek64.so.7.1"

// creation of an environment and task
MSKenv t env = NULL;
MSKtask t task = NULL;
MSKrescodee res;
/* Create an environment */
res = MSK makeenv(&env, NULL);
/* You may connect streams and other callbacks to env here */
/* Create a task */
if (res == MSK RES OK)
res = MSK maketask(env, 0,0, &task);
/* Load a problem into the task, optimize etc. */
MSK deletetask(&task);
MSK deleteenv(&env);

#define NUMCON 1668   
#define NUMVAR 2382  
#define NUMANZ 9236  
#define NUMQNZ 4764 


//gcc simple.c -o simple \ -I ~/mosek/7/tools/platform/linux64x86/h \  -L ~/mosek/7/tools/platform/linux64x86/bin \ -Wl,-rpath-link=~/mosek/7/tools/platform/linux64x86/platform/linux64x86/bin \ -lmosek64 -pthread