#pragma once
#include <stdio.h>

typedef struct file file;

struct file{
	int fd;
	fpos_t pos;
};

char *format_string(const char *format, ...);
void capture_out(char *line, char *cmd);
file toggle_capture_out_on(char *fpath);
void toggle_capture_out_off(file *f);
void load_file(char *path);


