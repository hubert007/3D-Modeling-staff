#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#ifdef __APPLE__
#include <GLEW/glew.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/freeglut.h>
#endif

float frand();
void read_shader_source_code(const char* filename, char *buf, int buf_size);

void check_shader_compile_error(GLuint shader);
void check_program_linking_error(GLuint program);
