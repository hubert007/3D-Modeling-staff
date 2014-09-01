#include "utils.h"
using namespace std;

float frand(){
	return (rand()%1000)/1000.0f;
}

void read_shader_source_code(const char* filename, char *buf, int buf_size)
{
    int file_size = 0;
    
    ifstream file(filename, ios::in);

    if( file ){
        file.seekg(0, ios::end);
        file_size = file.tellg();
        file.seekg (0, ios::beg);
    
        if(buf_size < file_size)
            cerr << "in read_shader_source_code: file size too big" << endl;
    
        file.read (buf, min(buf_size,file_size));
        file.close();
    
        buf[min(buf_size,file_size)] = '\0';

        file.close();
    }
    else{
        cerr << "can't open file " << filename << endl;
    }

    return;
    
}

void check_shader_compile_error(GLuint shader)
{
    GLint compile_status;
    char info[1024];
    int info_len = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &compile_status);

    if( !compile_status )
    {
        glGetShaderInfoLog(shader, 1024, &info_len, info);
        cerr << " Shader complie error: " << endl;
        cerr << info << endl;
    }
    else
    {
        cout << "compile shader successful" << endl;
    }
}

void check_program_linking_error(GLuint program)
{
    GLint link_status;
    char info[1024];
    int info_len = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &link_status);

    if( !link_status )
    {
        glGetProgramInfoLog(program, 1024, &info_len, info);
        cerr << " Program linking error: " << endl;
        cerr << info << endl;
    }
    else
    {
        cout << "link program successful" << endl;
    }
    
}

void screen_to_point(int x,int y,float d[2]){
	int viewport[4];
	double dX, dY, dClickY; // glUnProject uses doubles, but I'm using floats for these 3D vectors

	glGetIntegerv(GL_VIEWPORT, viewport);	
	int g_WindowHeight = glutGet(GLUT_WINDOW_HEIGHT);
	dClickY = double (g_WindowHeight - y); 

	d[0] = (float) dX;
	d[1] = (float) dY;
}
