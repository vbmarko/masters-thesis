#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>

 void capture_out(char *line,char * cmd)
    {
        int    fd;
        fpos_t pos;

        printf("stdout, ");

        fflush(stdout);
        fgetpos(stdout, &pos);
        fd = dup(fileno(stdout));
        freopen("stdout.out", "w", stdout);

        printf("%s",cmd);

        fflush(stdout);
        dup2(fd, fileno(stdout));
        close(fd);
        clearerr(stdout);
        fsetpos(stdout, &pos);        /* for C9X */
        printf("stdout again\n");
	
	FILE *fil = fopen("stdout.out","r");
	fgets(line,20,fil);
	fclose(fil);
    }

int main(){

	char* line = (char*) calloc(20,sizeof(char));
	capture_out(line, "this in line");
	printf("in the line we got %s",line);

	return 0;
}
