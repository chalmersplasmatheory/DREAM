#ifndef _DREAM_FVM_INIT_HPP
#define _DREAM_FVM_INIT_HPP

void dream_fvm_initialize();
void dream_fvm_finalize();

void dream_fvm_set_gsl_error_handler();
void dream_fvm_gsl_error_handler(const char*, const char*, int, int);

#endif/*_DREAM_FVM_INIT_HPP*/
