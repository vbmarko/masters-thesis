#include "nrn_lib/include/oc_ansi.h"
#include "nrn_lib/include/hocdec.h"
#include "nrn/src/nrniv/ocjump.h"
#include "nrn_lib/include/ocfunc.h"
//#include <bits/pthreadtypes.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "nrn_lib/include/hocdec.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <cstdarg>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>

#include "util.h"



extern int hoc_valid_stmt(const char* stmt, struct Object* ob);



char* format_string(const char* format, ...) {
    va_list args;
    va_start(args, format);

    // Calculate the required length for the formatted string
    int length = vsnprintf(NULL, 0, format, args);
    va_end(args);

    if (length < 0) {
        return NULL; // vsnprintf error
    }

    // Allocate memory for the formatted string
    char* result = (char*)malloc((length + 1) * sizeof(char));
    if (!result) {
        return NULL; // Memory allocation failed
    }

    // Format the string with the provided arguments
    va_start(args, format);
    vsnprintf(result, length + 1, format, args);
    va_end(args);

    return result;
}

 void capture_out(char *line,char * cmd)
    {
        int    fd;
        fpos_t pos;


        fflush(stdout);
        fgetpos(stdout, &pos);
        fd = dup(fileno(stdout));
        freopen("stdout.out", "w", stdout);

        hoc_valid_stmt(cmd,0);

        fflush(stdout);
        dup2(fd, fileno(stdout));
        close(fd);
        clearerr(stdout);
        fsetpos(stdout, &pos);        /* for C9X */

	
	FILE *fil = fopen("stdout.out","r");
	fgets(line,20,fil);
	fclose(fil);
    }

file toggle_capture_out_on(char *fpath)
    {
        int    fd;
        fpos_t pos;


        fflush(stdout);
        fgetpos(stdout, &pos);
        fd = dup(fileno(stdout));
	printf("from now print will go into the file %s. use with the returned file toggle_capture_out_off to print to stdoutagain",fpath);
        freopen(fpath, "w", stdout);
	return (file) {.fd = fd,.pos = pos};
    } 


void toggle_capture_out_off(file *f)
    {
        int    fd = f->fd;
        fpos_t pos = f->pos;


        fflush(stdout);
        dup2(fd, fileno(stdout));
        close(fd);
        clearerr(stdout);
        fsetpos(stdout, &pos);        /* for C9X */

	

    }

void load_file(char *path){
	char *cmd = format_string("load_file(\"%s\")", path);
	hoc_valid_stmt(cmd,0);
}
