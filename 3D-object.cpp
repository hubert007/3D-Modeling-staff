
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#ifdef __APPLE__
#include <GLEW/glew.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/freeglut.h>
#endif

#include "utils.h"
#include "debuglib.h"

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/transform.hpp"



// CONTSTANTS
const float PI = 3.1415926;
const int window_width = 500, window_height = 500;
const char * window_title = "3D-Object";

glm::vec3 POS;



struct Vertex
{
    glm::vec4 pos;
    glm::vec4 color;
    glm::vec4 normal;
};


// runtime variable
char buf[1024];
vector<Vertex> grid_vertices;


GLuint vertex_shader, fragment_shader, program;



// function Prototype 
void graphics_init();

void initialize_points(){
    // create the grid

    for(int i = 0; i < 11 ; i++)
    {
        Vertex start,end;
        start.color = glm::vec4(1,1,1,1);
        end.color = glm::vec4(1,1,1,1);

        start.pos = glm::vec4(i - 5.0, 0, -5.0, 1.0); 
        end.pos = glm::vec4(i - 5.0, 0, 5.0, 1.0);
        grid_vertices.push_back(start);
        grid_vertices.push_back(end);
            
        start.pos = glm::vec4(-5.0, 0, i  - 5.0, 1.0); 
        end.pos = glm::vec4(5.0f, 0, i - 5.0, 1.0); 
        grid_vertices.push_back(start);
        grid_vertices.push_back(end);
    }
    // create the frame at origin

    Vertex vo, vx,vy,vz;
    vo.pos = glm::vec4(0,0.01,0,1);
    
    vx.pos = glm::vec4(3,0.01,0,1);
    vx.color = glm::vec4(1,0,0,1); // red
    
    vy.pos = glm::vec4(0,3.01,0,1);
    vy.color = glm::vec4(0,1,0,1); // green
    
    vz.pos = glm::vec4(0,0.01,3,1);
    vz.color = glm::vec4(0,0,1,1); // green

    vo.color = vx.color;
    grid_vertices.push_back(vo);
    grid_vertices.push_back(vx);

    vo.color = vy.color;
    grid_vertices.push_back(vo);
    grid_vertices.push_back(vy);
    
    vo.color = vz.color;
    grid_vertices.push_back(vo);
    grid_vertices.push_back(vz); 
}

Vertex cuface[360];
Vertex cdface[360];
Vertex ccylinder[2];
Vertex ctriangle[360*3*2*2];
Vertex btriangle[360*3*2*2];
glm::mat4 bnewt= glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 0.25f, 0.0f));
void initialize_base(){
	//create base points at origin
	ccylinder[0].pos = glm::vec4(0,0.25,0,1);
	ccylinder[1].pos = glm::vec4(0,-0.25,0,1);
	
	for(int i=0; i<2; i++)
	{
	ccylinder[i].color = glm::vec4(0.5,0.5,0.5,1);
	}

	float dd = 0;
        
        //create upper and down part points of base
	for(int i=0; i<360; i++){
		cuface[i].pos = glm::vec4(sin(dd),0.25,cos(dd),1);
		cuface[i].color = glm::vec4(0.5,0.5,0.5,1);
		cdface[i].pos = glm::vec4(sin(dd),-0.25,cos(dd),1);
		cdface[i].color = glm::vec4(0.5,0.5,0.5,1);
		dd += 1.0f;
		}
	//normal	
	float x10 = cuface[0].pos.x - ccylinder[0].pos.x;
	float y10 = cuface[0].pos.y - ccylinder[0].pos.y;
	float z10 = cuface[0].pos.z - ccylinder[0].pos.z;
	float x11 = cuface[1].pos.x - ccylinder[0].pos.x;
	float y11 = cuface[1].pos.y - ccylinder[0].pos.y;
	float z11 = cuface[1].pos.z - ccylinder[0].pos.z;
	float x20 = cdface[0].pos.x - ccylinder[1].pos.x;
	float y20 = cdface[0].pos.y - ccylinder[1].pos.y;
	float z20 = cdface[0].pos.z - ccylinder[1].pos.z;
	float x21 = cdface[1].pos.x - ccylinder[1].pos.x;
	float y21 = cdface[1].pos.y - ccylinder[1].pos.y;
	float z21 = cdface[1].pos.z - ccylinder[1].pos.z;	
	ccylinder[0].normal = glm::normalize(glm::vec4((y10*z11-z10*y11),(z10*x11-x10*z11),(x10*y11-y10*x11),1.0f));
	ccylinder[1].normal = glm::normalize(glm::vec4((y21*z20-z21*y20),(z21*x20-x21*z20),(x21*y20-y21*x20),1.0f));

	glm::vec4 normalsPlane[360];
	for(int i=0; i<360; i++){
		int i_next = i+1;
		if(i==359) i_next=0;
		float xright = cuface[i_next].pos.x - cuface[i].pos.x;
		float yright = cuface[i_next].pos.y - cuface[i].pos.y;
		float zright = cuface[i_next].pos.z - cuface[i].pos.z;
		float xmiddle = cdface[i].pos.x - cuface[i].pos.x;
		float ymiddle = cdface[i].pos.y - cuface[i].pos.y;
		float zmiddle = cdface[i].pos.z - cuface[i].pos.z;
		glm::vec4 normal0 = glm::normalize(glm::vec4((ymiddle*zright-zmiddle*yright),(zmiddle*xright-xmiddle*zright),(xmiddle*yright-ymiddle*xright),1.0f));
		normalsPlane[i] = normal0;
	}	
	
	for(int i=0; i<360; i++){
		int i_before = i-1;
		if(i==0) i_before=359;
		cuface[i].normal = glm::normalize(ccylinder[0].normal + normalsPlane[i] + normalsPlane[i_before]);
		cdface[i].normal = glm::normalize(ccylinder[1].normal + normalsPlane[i] + normalsPlane[i_before]);
		}

	for(int i=0;i<360;i++){
		int i_next = i+1;
		if(i==359) 
                {
                 i_next=0;
                }
		ctriangle[i*3] = cuface[i];
		ctriangle[i*3+1] = ccylinder[0];
		ctriangle[i*3+2] = cuface[i_next];
		
		ctriangle[360*3+i*3] = cdface[i];
		ctriangle[360*3+i*3+1] = ccylinder[1];
		ctriangle[360*3+i*3+2] = cdface[i_next];
		
		ctriangle[360*6+i*3] = cuface[i];
		ctriangle[360*6+i*3+1] = cdface[i];
		ctriangle[360*6+i*3+2] = cuface[i_next];
		
		ctriangle[360*9+i*3] = cuface[i_next];
		ctriangle[360*9+i*3+1] = cdface[i];
		ctriangle[360*9+i*3+2] = cdface[i_next];
		}
        
        for(int i=0; i <360*3*2*2;i++)
	{
	btriangle[i].pos=bnewt*ctriangle[i].pos;
        btriangle[i].color= glm::vec4(0.5,0.5,0.5,1);
        btriangle[i].normal=ctriangle[i].normal;
	}
	
	}

glm::mat4 tnews=glm::scale(glm::mat4(1.0f),glm::vec3(0.8f,1.6f,0.8f));
glm::mat4 tnewt= glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 1.25f, 0.0f));
glm::mat4 tnewr = glm::rotate(glm::mat4(1.0f),90.0f,glm::vec3(1.0f,0.0f,0.0f));
Vertex ttriangle[360*3*2*2];
void initialize_top(){
	for(int i=0; i <360*3*2*2;i++)
	{	
        ttriangle[i].color = glm::vec4(0,0,1,1);
	}
	
        for(int i=0; i <360*3*2*2;i++)
	{
	ttriangle[i].pos=tnewt*tnewr*tnews*ctriangle[i].pos;	
        ttriangle[i].normal=tnewr*ctriangle[i].normal;	
	}
	}

Vertex cctriangle[6*2*3];
Vertex cube[8];
Vertex ftriangle[6*2*3];
glm::mat4 fnews=glm::scale(glm::mat4(1.0f),glm::vec3(0.5f,0.2f,0.1f));
glm::mat4 fnewt= glm::translate(glm::mat4(1.0f),glm::vec3(1.1f, 1.2f, 0.0f));
glm::mat4 fnewt1=glm::translate(glm::mat4(1.0f),glm::vec3(1.1f, 0.0f, 0.0f));
glm::mat4 fnewt2=glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 1.2f, 0.0f));
void initialize_firstarm(){
    //create cube points at origin
    cube[0].pos = glm::vec4(1,1,1,1);
    cube[1].pos = glm::vec4(1,1,-1,1);
    cube[2].pos = glm::vec4(1,-1,1,1);
    cube[3].pos = glm::vec4(1,-1,-1,1);
    cube[4].pos = glm::vec4(-1,1,1,1);
    cube[5].pos = glm::vec4(-1,1,-1,1);
    cube[6].pos = glm::vec4(-1,-1,1,1);
    cube[7].pos = glm::vec4(-1,-1,-1,1);
    
     for(int i =0;i<8;i++){
	cube[i].color = glm::vec4(1,1,0,1);
	}

//calculate normal
//normal	
	float x10 = cube[1].pos.x - cube[0].pos.x; 
	float y10 = cube[1].pos.y - cube[0].pos.y;
	float z10 = cube[1].pos.z - cube[0].pos.z;
	float x20 = cube[2].pos.x - cube[0].pos.x; 
	float y20 = cube[2].pos.y - cube[0].pos.y;
	float z20 = cube[2].pos.z - cube[0].pos.z;
	float x40 = cube[4].pos.x - cube[0].pos.x; 
	float y40 = cube[4].pos.y - cube[0].pos.y;
	float z40 = cube[4].pos.z - cube[0].pos.z;
	float x37 = cube[3].pos.x - cube[7].pos.x; 
	float y37 = cube[3].pos.y - cube[7].pos.y;
	float z37 = cube[3].pos.z - cube[7].pos.z;
	float x57 = cube[5].pos.x - cube[7].pos.x; 
	float y57 = cube[5].pos.y - cube[7].pos.y;
	float z57 = cube[5].pos.z - cube[7].pos.z;
	float x67 = cube[6].pos.x - cube[7].pos.x; 
	float y67 = cube[6].pos.y - cube[7].pos.y;
	float z67 = cube[6].pos.z - cube[7].pos.z;	
	glm::vec3 normal012 = glm::normalize(glm::vec3((y20*z10-z20*y10),(z20*x10-x20*z10),(x20*y10-y20*x10)));//plane 0123
	glm::vec3 normal014 = glm::normalize(glm::vec3((y10*z40-z10*y40),(z10*x40-x10*z40),(x10*y40-y10*x40)));//plane 0145
	glm::vec3 normal024 = glm::normalize(glm::vec3((y40*z20-z40*y20),(z40*x20-x40*z20),(x40*y20-y40*x20)));//plane 0246
	glm::vec3 normal735 = glm::normalize(glm::vec3((y57*z37-z57*y37),(z57*x37-x57*z37),(x57*y37-y57*x37)));//plane 1357
	glm::vec3 normal736 = glm::normalize(glm::vec3((y37*z67-z37*y67),(z37*x67-x37*z67),(x37*y67-y37*x67)));//plane 2367	
	glm::vec3 normal756 = glm::normalize(glm::vec3((y67*z57-z67*y57),(z67*x57-x67*z57),(x67*y57-y67*x57)));//plane 4567
	
	cube[0].normal = glm::vec4(glm::normalize(normal012 + normal014 + normal024),1.0f);
	cube[1].normal = glm::vec4(glm::normalize(normal012 + normal014 + normal735),1.0f);
	cube[2].normal = glm::vec4(glm::normalize(normal012 + normal024 + normal736),1.0f);
	cube[3].normal = glm::vec4(glm::normalize(normal012 + normal735 + normal736),1.0f);
	cube[4].normal = glm::vec4(glm::normalize(normal756 + normal014 + normal024),1.0f);
	cube[5].normal = glm::vec4(glm::normalize(normal756 + normal014 + normal735),1.0f);
	cube[6].normal = glm::vec4(glm::normalize(normal756 + normal024 + normal736),1.0f);
	cube[7].normal = glm::vec4(glm::normalize(normal756 + normal736 + normal735),1.0f);	 




    //the pos of 12 triangles	    
    cctriangle[0] =cctriangle[3] =cctriangle[6] =cctriangle[9] =cctriangle[12] =cctriangle[15] = cube[0];
    cctriangle[1] =cctriangle[7] =cctriangle[20] =cctriangle[23] = cube[1];
    cctriangle[4] =cctriangle[16] =cctriangle[32] =cctriangle[35] =cube[2];
    cctriangle[2] =cctriangle[5] =cctriangle[22] =cctriangle[31] =cube[3];
    cctriangle[10] =cctriangle[13] =cctriangle[26] =cctriangle[29] =cube[4];
    cctriangle[8] =cctriangle[11] =cctriangle[19] =cctriangle[25] =cube[5];
    cctriangle[14] =cctriangle[17] =cctriangle[28] =cctriangle[34] =cube[6];
    cctriangle[18] =cctriangle[21] =cctriangle[24] =cctriangle[27] =cctriangle[30] =cctriangle[33] = cube[7];

   
    for(int i=0; i <36;i++)
	{
        ftriangle[i].color=glm::vec4(1,1,0,1);
	ftriangle[i].pos=fnewt*fnews*cctriangle[i].pos;	
        ftriangle[i].normal=cctriangle[i].normal;
	}
}
	
glm::mat4 jnews=glm::scale(glm::mat4(1.0f),glm::vec3(0.3f,0.8f,0.3f));
glm::mat4 jnewr = glm::rotate(glm::mat4(1.0f),90.0f,glm::vec3(1.0f,0.0f,0.0f));
glm::mat4 jnewt= glm::translate(glm::mat4(1.0f),glm::vec3(1.7f, 1.2f, 0.0f));  
glm::mat4 jnewt1= glm::translate(glm::mat4(1.0f),glm::vec3(1.7f, 0.0f, 0.0f)); 
glm::mat4 jnewt2= glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 1.2f, 0.0f)); 
Vertex jtriangle[360*3*2*2];
void initialize_joint(){
	for(int i=0; i <360*3*2*2;i++)
	{	
        jtriangle[i].color = glm::vec4(0,1,1,1);
	}
	for(int i=0; i <360*3*2*2;i++)
	{
	jtriangle[i].pos=jnewt*jnewr*jnews*ctriangle[i].pos;
        jtriangle[i].normal=jnewr*ctriangle[i].normal;	
	}
	}


glm::mat4 snews=glm::scale(glm::mat4(1.0f),glm::vec3(0.05f,1.8f,0.05f));
glm::mat4 snewr = glm::rotate(glm::mat4(1.0f),90.0f,glm::vec3(0.0f,0.0f,1.0f));
glm::mat4 snewt= glm::translate(glm::mat4(1.0f),glm::vec3(2.25f, 1.2f, 0.0f));
glm::mat4 snewt1= glm::translate(glm::mat4(1.0f),glm::vec3(2.25f, 0.0f, 0.0f));
glm::mat4 sjnewt1= glm::translate(glm::mat4(1.0f),glm::vec3(0.55f, 0.0f, 0.0f));
glm::mat4 snewt2= glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 1.2f, 0.0f));
Vertex striangle[360*3*2*2];
void initialize_secondarm(){
	for(int i=0; i <360*3*2*2;i++)
	{	
        striangle[i].color = glm::vec4(1,0,1,1);
	}
	
        for(int i=0; i <360*3*2*2;i++)
	{
	striangle[i].pos=snewt*snewr*snews*ctriangle[i].pos;	
        striangle[i].normal=snewr*ctriangle[i].pos;
	}
	}


glm::mat4 pnews=glm::scale(glm::mat4(1.0f),glm::vec3(0.05f,1.8f,0.05f));
glm::mat4 pnewr = glm::rotate(glm::mat4(1.0f),90.0f,glm::vec3(0.0f,0.0f,1.0f));
glm::mat4 pnewt= glm::translate(glm::mat4(1.0f),glm::vec3(3.15f, 1.2f, 0.0f));
glm::mat4 pnewt1= glm::translate(glm::mat4(1.0f),glm::vec3(3.15f, 0.0f, 0.0f));
glm::mat4 pxnewt1= glm::translate(glm::mat4(1.0f),glm::vec3(0.45f, 0.0f, 0.0f));
glm::mat4 pxnewt2= glm::translate(glm::mat4(1.0f),glm::vec3(1.0f, 0.0f, 0.0f));
glm::mat4 pjnewt1= glm::translate(glm::mat4(1.0f),glm::vec3(1.45f, 0.0f, 0.0f));
glm::mat4 psnewt1= glm::translate(glm::mat4(1.0f),glm::vec3(0.9f, 0.0f, 0.0f));
glm::mat4 pnewt2= glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 1.2f, 0.0f));
Vertex ptriangle[360*3*2*2];
void initialize_pen(){
	for(int i=0; i <360*3*2*2;i++)
	{	
        ptriangle[i].color = glm::vec4(0.5,1,0.5,1);
	}
	
        for(int i=0; i <360*3*2*2;i++)
	{
	ptriangle[i].pos=pnewt*pnewr*pnews*ctriangle[i].pos;
        ptriangle[i].normal=pnewr*ctriangle[i].normal;	
	}
	}


glm::mat4 xnews=glm::scale(glm::mat4(1.0f),glm::vec3(0.05f));
glm::mat4 xnewt= glm::translate(glm::mat4(1.0f),glm::vec3(3.4f, 1.25f, 0.0f));
glm::mat4 xnewt1= glm::translate(glm::mat4(1.0f),glm::vec3(3.4f, 0.0f, 0.0f));
glm::mat4 xjnewt1=glm::translate(glm::mat4(1.0f),glm::vec3(1.7f, 0.0f, 0.0f));
glm::mat4 xsnewt1=glm::translate(glm::mat4(1.0f),glm::vec3(0.7f, 0.0f, 0.0f));
glm::mat4 xsnewt2=glm::translate(glm::mat4(1.0f),glm::vec3(1.0f, 0.0f, 0.0f));
glm::mat4 xnewt2= glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 1.25f, 0.0f));
Vertex xtriangle[6*2*3];
void initialize_box(){
	    for(int i=0; i <36;i++)
	{
        xtriangle[i].color=glm::vec4(1,1,0,1);
	xtriangle[i].pos=xnewt*xnews*cctriangle[i].pos;
        xtriangle[i].normal=cctriangle[i].normal;	
	}
}

glm::mat4 hnews=glm::scale(glm::mat4(1.0f),glm::vec3(0.25f,0.6f,0.5f));
glm::mat4 hnewt= glm::translate(glm::mat4(1.0f),glm::vec3(-2.0f, 3.0f, 0.0f));
Vertex htriangle[6*2*3];
void initialize_head()
{
 for(int i=0; i <36;i++)
	{
        htriangle[i].color=glm::vec4(1,1,0.5,1);
	htriangle[i].pos=hnewt*hnews*cctriangle[i].pos;
        htriangle[i].normal=cctriangle[i].normal;	
	}
}

glm::mat4 bonews=glm::scale(glm::mat4(1.0f),glm::vec3(0.3f,1.2f,0.8f));
glm::mat4 bonewt= glm::translate(glm::mat4(1.0f),glm::vec3(-2.0f,1.7f, 0.0f));
Vertex botriangle[6*2*3];
void initialize_body()
{
for(int i=0; i <36;i++)
	{
        botriangle[i].color=glm::vec4(1,0.2,1,1);
	botriangle[i].pos=bonewt*bonews*cctriangle[i].pos;
        botriangle[i].normal=cctriangle[i].normal;	
	}
}

glm::mat4 anews1=glm::scale(glm::mat4(1.0f),glm::vec3(0.8f,0.2f,0.2f));
glm::mat4 anewt1= glm::translate(glm::mat4(1.0f),glm::vec3(-1.0f,1.5f, 0.8f));
glm::mat4 anewr1 = glm::rotate(anewr1,-10.0f,glm::vec3(0.0f,0.0f,1.0f));  
Vertex atriangle1[6*2*3];
void initialize_arm1()
{
for(int i=0; i <36;i++)
	{
        atriangle1[i].color=glm::vec4(1,1,0.2,1);
	atriangle1[i].pos=anewt1*anews1*cctriangle[i].pos;
        atriangle1[i].normal=cctriangle[i].normal;	
	}
}


glm::mat4 anews2=glm::scale(glm::mat4(1.0f),glm::vec3(0.8f,0.2f,0.2f));
glm::mat4 anewt2= glm::translate(glm::mat4(1.0f),glm::vec3(-1.0f,1.5f, -0.8f));
glm::mat4 anewr2 = glm::rotate(anewr1,-45.0f,glm::vec3(0.0f,0.0f,1.0f));  
Vertex atriangle2[6*2*3];
void initialize_arm2()
{
for(int i=0; i <36;i++)
	{
        atriangle2[i].color=glm::vec4(1,1,0.2,1);
	atriangle2[i].pos=anewt2*anews2*cctriangle[i].pos;
        atriangle2[i].normal=cctriangle[i].normal;	
	}
}

glm::mat4 lnews1=glm::scale(glm::mat4(1.0f),glm::vec3(0.2f,1.2f,0.25f));
glm::mat4 lnewt1= glm::translate(glm::mat4(1.0f),glm::vec3(-2.0f,0.6f, 0.4f));
Vertex ltriangle1[6*2*3];
void initialize_leg1()
{
for(int i=0; i <36;i++)
	{
        ltriangle1[i].color=glm::vec4(1,1,0.2,1);
	ltriangle1[i].pos=lnewt1*lnews1*cctriangle[i].pos;
        ltriangle1[i].normal=cctriangle[i].normal;	
	}
}

glm::mat4 lnews2=glm::scale(glm::mat4(1.0f),glm::vec3(0.2f,1.2f,0.25f));
glm::mat4 lnewt2= glm::translate(glm::mat4(1.0f),glm::vec3(-2.0f,0.6f, -0.4f));
Vertex ltriangle2[6*2*3];
void initialize_leg2()
{
for(int i=0; i <36;i++)
	{
        ltriangle2[i].color=glm::vec4(1,1,0.2,1);
	ltriangle2[i].pos=lnewt2*lnews2*cctriangle[i].pos;
        ltriangle2[i].normal=cctriangle[i].normal;	
	}
}

// MOUSE handling *******************************************
int last_x, last_y;
int selected_idx = -1;

void mouse(int button, int state, int x, int y ){
	if(state == GLUT_DOWN){
        
	}else{
        
	}

    last_x = x;
    last_y = y;
    
}


void motion( int x, int y){
    GLint viewport[4];

    glutPostRedisplay();
    
}

int hint=0;
int base=0;
int top=0;
int first=0;
int second=0;
int pen=0;
int light=0;
int bonus=0;
// KEYBOARD handling *******************************************
void keyboard(unsigned char key, int x, int y)
{
    switch (key) {
    case 'r':
        graphics_init();
        glutPostRedisplay();
        break;
    case 27:
        exit(0);
        break; 
    case 99:
        if(hint==0)
        {  hint=1;
          cout<<"camera is selected"<<endl;
        }
        else if (hint=1)
        {  hint=0; 
          cout<<"camera is not selected"<<endl;
        }
     break;
    case 49:
	if(base==0)
        {  base=1;
           hint==0;
          cout<<"model is selected"<<endl;
        }
        else if (base==1)
        {  base=0; 
          cout<<"model is not selected"<<endl;
        }
     break;
     case 50:
	if(top==0)
        {  top=1;
           hint=0;
           base=0;
          cout<<"top is selected"<<endl;
        }
        else if (top==1)
        {  top=0; 
          cout<<"top is not selected"<<endl;
        }
     break;
     case 51:
	if(first==0)
        {  first=1;
           hint=0;
           base=0;
           top=0;
          cout<<"first arm is selected"<<endl;
        }
        else if (first==1)
        {  first=0; 
          cout<<"first arm is not selected"<<endl;
        }
     break;
     case 52:
	if(second==0)
        {  second=1;
           hint=0;
           base=0;
           top=0;
           first=0;
          cout<<"second arm is selected"<<endl;
        }
        else if (second==1)
        {  second=0; 
          cout<<"second arm is not selected"<<endl;
        }
     break;
     case 54:
        if(bonus==0)
        {  bonus=1;
           cout<<"bonus is selected"<<endl;
        }
        else if(bonus==1)
        {  bonus=0;
           cout<<"bonus is not selected"<<endl;
        }
     break;
     case 53:
	if(pen==0)
        {  pen=1;
           hint=0;
           base=0;
           top=0;
           second=0;
          cout<<"pen is selected"<<endl;
        }
        else if (pen==1)
        {  pen=0; 
          cout<<"pen is not selected"<<endl;
        }
     break;
     
     case 'l':
        if(light==0)
        {  light=1;
           pen=0;
           hint=0;
           base=0;
           top=0;
           second=0;
          cout<<"light is on"<<endl;
        }
        else if (light==1)
        {  light=0; 
           cout<<"light is off"<<endl;
        }
      
    glutPostRedisplay();}
}

glm::mat4 new1;
glm::mat4 new2;
glm::mat4 new3;
glm::mat4 new4;
glm::mat4 new5;
glm::mat4 new6;
glm::mat4 newt; 
float cx=0;float cy=0;
void special_key(int key, int x, int y)
{

    // TODO: capture the arrow key here  
   if(hint==1)
	{ if(key==GLUT_KEY_LEFT)
          {
              cy=cy+0.05f;
           
	  }
 	  else if(key==GLUT_KEY_RIGHT)
          {
	      cy=cy-0.05f;
	  }
 	  else if(key==GLUT_KEY_UP)
          {
              cx=cx+0.05f;
	  }
          else if(key==GLUT_KEY_DOWN)
          {
 	      cx=cx-0.05f;
	  }
        }
   else if(base == 1){
     
     if(key==GLUT_KEY_LEFT)
          {
          for(int i=0;i<36;i++){
          new1 = glm::translate(glm::mat4(1.0f),glm::vec3(-1.0f, 0.0f, 0.0f));
	  ftriangle[i].pos = new1 *ftriangle[i].pos;
          xtriangle[i].pos = new1 *xtriangle[i].pos;   
          htriangle[i].pos = new1 *htriangle[i].pos;
	  botriangle[i].pos = new1 *botriangle[i].pos;
	  atriangle1[i].pos = new1 *atriangle1[i].pos;
	  atriangle2[i].pos = new1 *atriangle2[i].pos;
	  ltriangle1[i].pos = new1 *ltriangle1[i].pos;     
	  ltriangle2[i].pos = new1 *ltriangle2[i].pos;	  
	}
	  for(int i=0;i<360*12;i++){
	  new1 = glm::translate(glm::mat4(1.0f),glm::vec3(-1.0f, 0.0f, 0.0f));
	  btriangle[i].pos = new1 * btriangle[i].pos;
          ttriangle[i].pos = new1 * ttriangle[i].pos;
          jtriangle[i].pos = new1 * jtriangle[i].pos;
          striangle[i].pos = new1 * striangle[i].pos;
          ptriangle[i].pos = new1 * ptriangle[i].pos;            
	  }
          }
     else if(key==GLUT_KEY_RIGHT)
          {
	  for(int i=0;i<36;i++){
          new1 = glm::translate(glm::mat4(1.0f),glm::vec3(1.0f, 0.0f, 0.0f));
	  ftriangle[i].pos = new1 *ftriangle[i].pos;
	  xtriangle[i].pos = new1 *xtriangle[i].pos;
	  htriangle[i].pos = new1 *htriangle[i].pos;
	  botriangle[i].pos = new1 *botriangle[i].pos;
	  atriangle1[i].pos = new1 *atriangle1[i].pos;
	  atriangle2[i].pos = new1 *atriangle2[i].pos;
	  ltriangle1[i].pos = new1 *ltriangle1[i].pos;     
	  ltriangle2[i].pos = new1 *ltriangle2[i].pos;
	  }
          for(int i=0;i<360*12;i++){
	  new1 = glm::translate(glm::mat4(1.0f),glm::vec3(1.0f, 0.0f, 0.0f));
	  btriangle[i].pos = new1 * btriangle[i].pos;
          ttriangle[i].pos = new1 * ttriangle[i].pos;
          jtriangle[i].pos = new1 * jtriangle[i].pos;
          striangle[i].pos = new1 * striangle[i].pos;
          ptriangle[i].pos = new1 * ptriangle[i].pos;            
	  }
 	  }
     else if(key==GLUT_KEY_UP)
          {
          for(int i=0;i<36;i++){
          new1 = glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 0.0f, -1.0f));
	  ftriangle[i].pos = new1 *ftriangle[i].pos;
	  xtriangle[i].pos = new1 *xtriangle[i].pos;
	  htriangle[i].pos = new1 *htriangle[i].pos;
	  botriangle[i].pos = new1 *botriangle[i].pos;
	  atriangle1[i].pos = new1 *atriangle1[i].pos;
	  atriangle2[i].pos = new1 *atriangle2[i].pos;
	  ltriangle1[i].pos = new1 *ltriangle1[i].pos;     
	  ltriangle2[i].pos = new1 *ltriangle2[i].pos;
	  }
	  for(int i=0;i<360*12;i++){
	  new1 = glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 0.0f, -1.0f));
	  btriangle[i].pos = new1 * btriangle[i].pos;
          ttriangle[i].pos = new1 * ttriangle[i].pos;
          jtriangle[i].pos = new1 * jtriangle[i].pos;
          striangle[i].pos = new1 * striangle[i].pos;
          ptriangle[i].pos = new1 * ptriangle[i].pos;            
	  }
          }
     else if(key==GLUT_KEY_DOWN)
          {
 	  for(int i=0;i<36;i++){
          new1 = glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 0.0f, 1.0f));
	  ftriangle[i].pos = new1 *ftriangle[i].pos;
          xtriangle[i].pos = new1 *xtriangle[i].pos;
          htriangle[i].pos = new1 *htriangle[i].pos;
	  botriangle[i].pos = new1 *botriangle[i].pos;
	  atriangle1[i].pos = new1 *atriangle1[i].pos;
	  atriangle2[i].pos = new1 *atriangle2[i].pos;
	  ltriangle1[i].pos = new1 *ltriangle1[i].pos;     
	  ltriangle2[i].pos = new1 *ltriangle2[i].pos;
	  }
	  for(int i=0;i<360*12;i++){
	  new1 = glm::translate(glm::mat4(1.0f),glm::vec3(0.0f, 0.0f, 1.0f));
	  btriangle[i].pos = new1 * btriangle[i].pos;
          ttriangle[i].pos = new1 * ttriangle[i].pos;
          jtriangle[i].pos = new1 * jtriangle[i].pos;
          striangle[i].pos = new1 * striangle[i].pos;
          ptriangle[i].pos = new1 * ptriangle[i].pos;            
	  }
	  }
     }
     else if(top==1)
    {
          if(key==GLUT_KEY_LEFT) 
          {
	  new2 = glm::rotate(new2,20.0f,glm::vec3(0.0f,1.0f,0.0f));   
	  }
	  else if(key==GLUT_KEY_RIGHT)
	  {
	  new2 = glm::rotate(new2,-20.0f,glm::vec3(0.0f,1.0f,0.0f));
	  }
	  float xd=btriangle[1].pos[0];
          float zd=btriangle[1].pos[2];
          newt = glm::translate(glm::mat4(1.0f),glm::vec3(xd,0.0f,zd));
          for(int i=0;i<36;i++){
	  ftriangle[i].pos = newt*new2*fnewt*fnews*cctriangle[i].pos;
          ftriangle[i].normal=new2*fnewt*cctriangle[i].normal;
          xtriangle[i].pos = newt*new2*xnewt*xnews*cctriangle[i].pos;  
	  xtriangle[i].normal=new2*xnewt*cctriangle[i].normal;        
	  }
	  for(int i=0;i<360*12;i++){	     
          
          ttriangle[i].pos =   newt*new2*tnewt*tnewr*tnews* ctriangle[i].pos;
          ttriangle[i].normal=new2*tnewt*tnewr*ctriangle[i].normal;
          jtriangle[i].pos =   newt*new2*jnewt*jnewr*jnews* ctriangle[i].pos;
	  jtriangle[i].normal=new2*jnewt*jnewr*ctriangle[i].normal;
          striangle[i].pos =   newt*new2*snewt*snewr*snews* ctriangle[i].pos;
          striangle[i].normal=new2*snewt*snewr*ctriangle[i].normal;
          ptriangle[i].pos =   newt*new2*pnewt*pnewr*pnews* ctriangle[i].pos;            
	  ptriangle[i].normal=new2*pnewt*pnewr*ctriangle[i].normal;
	  }
    }
    else if(first==1)
    {
	  if(key==GLUT_KEY_UP) 
          {
	  new3 = glm::rotate(new3,20.0f,glm::vec3(0.0f,0.0f,1.0f));   
	  }
	  else if(key==GLUT_KEY_DOWN)
	  {
	  new3 = glm::rotate(new3,-20.0f,glm::vec3(0.0f,0.0f,1.0f));
	  }
	  float xd=btriangle[1].pos[0];
          float zd=btriangle[1].pos[2];
          newt = glm::translate(glm::mat4(1.0f),glm::vec3(xd,0.0f,zd));
	  for(int i=0;i<36;i++){
	  ftriangle[i].pos = newt*fnewt2*new2*new3*fnewt1*fnews*cctriangle[i].pos;
          ftriangle[i].normal=new2*new3*fnewt1*cctriangle[i].normal;
          xtriangle[i].pos = newt*xnewt2*new2*new3*xnewt1*xnews*cctriangle[i].pos; 
          xtriangle[i].normal=new2*new3*xnewt1*cctriangle[i].normal;
	  }
	  for(int i=0;i<360*12;i++){	     
          jtriangle[i].pos =   newt*jnewt2*new2*new3*jnewt1*jnewr*jnews* ctriangle[i].pos;
          jtriangle[i].normal=new2*new3*jnewt1*jnewr*ctriangle[i].normal;
          striangle[i].pos =   newt*snewt2*new2*new3*snewt1*snewr*snews* ctriangle[i].pos;
          striangle[i].normal=new2*new3*snewt1*snewr*ctriangle[i].normal;
          ptriangle[i].pos =   newt*pnewt2*new2*new3*pnewt1*pnewr*pnews* ctriangle[i].pos;   
	  ptriangle[i].normal=new2*new3*pnewt1*pnewr*ctriangle[i].normal;         
	  }
      
    }
    else if(second==1)
    {
	if(key==GLUT_KEY_UP) 
          {
	  new4 = glm::rotate(new4,20.0f,glm::vec3(0.0f,0.0f,1.0f));   
	  }
	  else if(key==GLUT_KEY_DOWN)
	  {
	  new4 = glm::rotate(new4,-20.0f,glm::vec3(0.0f,0.0f,1.0f));
	  }
	  float xd=btriangle[1].pos[0];
          float zd=btriangle[1].pos[2];
          newt = glm::translate(glm::mat4(1.0f),glm::vec3(xd,0.0f,zd));
	  for(int i=0;i<36;i++){
          xtriangle[i].pos = newt*xnewt2*new2*new3*jnewt1*new4*xjnewt1*xnews*cctriangle[i].pos;  
	  xtriangle[i].normal=new2*new3*jnewt1*new4*xjnewt1*cctriangle[i].normal;        
	  }
	  for(int i=0;i<360*12;i++){	     
          striangle[i].pos =   newt*snewt2*new2*new3*jnewt1*new4*sjnewt1*snewr*snews* ctriangle[i].pos;
	  striangle[i].normal=new2*new3*jnewt1*new4*sjnewt1*snewr*ctriangle[i].normal;
          ptriangle[i].pos =   newt*pnewt2*new2*new3*jnewt1*new4*pjnewt1*pnewr*pnews* ctriangle[i].pos;
	  ptriangle[i].normal=new2*new3*jnewt1*new4*pjnewt1*pnewr*ctriangle[i].normal;           
	  }
    }
    else if(pen==1)
    {
     float xd=btriangle[1].pos[0];
     float zd=btriangle[1].pos[2];
     newt = glm::translate(glm::mat4(1.0f),glm::vec3(xd,0.0f,zd));
     if(key==GLUT_KEY_UP) 
          {
	  new5 = glm::rotate(new5,20.0f,glm::vec3(0.0f,0.0f,1.0f));
	  for(int i=0;i<36;i++){
          xtriangle[i].pos = newt*xnewt2*new2*new3*jnewt1*new4*xsnewt2*new5*xsnewt1*xnews*cctriangle[i].pos;       
	  xtriangle[i].normal=new2*new3*jnewt1*new4*xsnewt2*new5*xsnewt1*cctriangle[i].normal;   
	  }
	  for(int i=0;i<360*12;i++){	     
          ptriangle[i].pos =   newt*pnewt2*new2*new3*jnewt1*new4*pxnewt2*new5*pxnewt1*pnewr*pnews* ctriangle[i].pos;
	  ptriangle[i].normal=new2*new3*jnewt1*new4*pxnewt2*new5*pxnewt1*pnewr* ctriangle[i].normal;
	  }   
	  }
     else if(key==GLUT_KEY_DOWN)
	  {
	  new5 = glm::rotate(new5,-20.0f,glm::vec3(0.0f,0.0f,1.0f));
	  for(int i=0;i<36;i++){
          xtriangle[i].pos = newt*xnewt2*new2*new3*jnewt1*new4*xsnewt2*new5*xsnewt1*xnews*cctriangle[i].pos;        
	  xtriangle[i].normal=new2*new3*jnewt1*new4*xsnewt2*new5*xsnewt1*cctriangle[i].normal;   
	  }
	  for(int i=0;i<360*12;i++){	     
          ptriangle[i].pos =   newt*pnewt2*new2*new3*jnewt1*new4*pxnewt2*new5*pxnewt1*pnewr*pnews* ctriangle[i].pos;
	  ptriangle[i].normal=new2*new3*jnewt1*new4*pxnewt2*new5*pxnewt1*pnewr* ctriangle[i].normal;            
	  }
	  }
     else if(glutGetModifiers()==GLUT_ACTIVE_SHIFT && key==GLUT_KEY_LEFT)
	{
        new6 = glm::rotate(new6,-20.0f,glm::vec3(1.0f,0.0f,0.0f));
	  for(int i=0;i<36;i++){
          xtriangle[i].pos = newt*xnewt2*new2*new3*jnewt1*new4*xsnewt2*new6*new5*xsnewt1*xnews*cctriangle[i].pos;  
	  xtriangle[i].normal=new2*new3*jnewt1*new4*xsnewt2*new6*new5*xsnewt1*cctriangle[i].normal;        
	  }
	  for(int i=0;i<360*12;i++){	     
          ptriangle[i].pos =   newt*pnewt2*new2*new3*jnewt1*new4*pxnewt2*new6*new5*pxnewt1*pnewr*pnews* ctriangle[i].pos;
	  ptriangle[i].normal=new2*new3*jnewt1*new4*pxnewt2*new6*new5*pxnewt1*pnewr* ctriangle[i].normal;            
	  }
	}
     else if(glutGetModifiers()==GLUT_ACTIVE_SHIFT && key==GLUT_KEY_RIGHT)
	{
       new6 = glm::rotate(new6,20.0f,glm::vec3(1.0f,0.0f,0.0f));
	  for(int i=0;i<36;i++){
          xtriangle[i].pos = newt*xnewt2*new2*new3*jnewt1*new4*xsnewt2*new6*new5*xsnewt1*xnews*cctriangle[i].pos;   
	  xtriangle[i].normal=new2*new3*jnewt1*new4*xsnewt2*new6*new5*xsnewt1*cctriangle[i].normal;         
	  }
	  for(int i=0;i<360*12;i++){	     
          ptriangle[i].pos =   newt*pnewt2*new2*new3*jnewt1*new4*pxnewt2*new6*new5*pxnewt1*pnewr*pnews* ctriangle[i].pos;            
	  }
	}
     else if(key==GLUT_KEY_LEFT)
	  {
	  new5 = glm::rotate(new5,-20.0f,glm::vec3(0.0f,1.0f,0.0f));
          for(int i=0;i<36;i++){
          xtriangle[i].pos = newt*xnewt2*new2*new3*jnewt1*new4*xsnewt2*new5*xsnewt1*xnews*cctriangle[i].pos;  
	  xtriangle[i].normal=new2*new3*jnewt1*new4*xsnewt2*new5*xsnewt1*cctriangle[i].normal;         
	  }
	  for(int i=0;i<360*12;i++){	     
          ptriangle[i].pos =   newt*pnewt2*new2*new3*jnewt1*new4*pxnewt2*new5*pxnewt1*pnewr*pnews* ctriangle[i].pos;            
	  ptriangle[i].normal=new2*new3*jnewt1*new4*pxnewt2*new5*pxnewt1*pnewr* ctriangle[i].normal;	  
	}
	  }
     else if(key==GLUT_KEY_RIGHT)
	  {
	  new5 = glm::rotate(new5,20.0f,glm::vec3(0.0f,1.0f,0.0f));
	for(int i=0;i<36;i++){
          xtriangle[i].pos = newt*xnewt2*new2*new3*jnewt1*new4*xsnewt2*new5*xsnewt1*xnews*cctriangle[i].pos;  
	  xtriangle[i].normal=new2*new3*jnewt1*new4*xsnewt2*new5*xsnewt1*cctriangle[i].normal;         
	  }
	for(int i=0;i<360*12;i++){	     
         ptriangle[i].pos =   newt*pnewt2*new2*new3*jnewt1*new4*pxnewt2*new5*pxnewt1*pnewr*pnews* ctriangle[i].pos;            
	 ptriangle[i].normal=new2*new3*jnewt1*new4*pxnewt2*new5*pxnewt1*pnewr* ctriangle[i].normal; 
	 }   
	  }
    }
    

    glutPostRedisplay();

}

// DISPLAY and RENDERING functions *************************

void lighting()
{
//intensity=emission+(1/(0+0*d+0.5*d*d))*spotl((p-v)/d*(1,1,1,1)diffusem+(s*n)^shiness*(1,1,1,1)*specularm

	glm::vec4 light = glm::vec4(10.0f,10.0f,10.0f,1.0f);
	glm::vec4 diffuse = glm::vec4(1,1,1,1);
        glm::vec4 diffusem=glm::vec4(0.8,0.8,0.8,1);
	glm::vec4 specular = glm::vec4(1,1,1,1);
        glm::vec4 specularm = glm::vec4(0.8,0.8,0.8,1);
        float diffusex=diffuse.x*diffusem.x;
        float diffusey=diffuse.y*diffusem.y;
        float diffusez=diffuse.z*diffusem.z;
        float specularx=specular.x*specularm.x;
        float speculary=specular.y*specularm.y;
        float specularz=specular.z*specularm.z;
//calculate base
	for(int i = 0; i<360*12; i++){
/*        float d2 = ((light.x-btriangle[i].pos.x)*(light.x-btriangle[i].pos.x)+(light.y-btriangle[i].pos.y)*(light.y-btriangle[i].pos.y)+(light.z-btriangle[i].pos.z)*(light.z-btriangle[i].pos.z));*/
	float diffuse_intensity = max(glm::dot(glm::normalize(light - btriangle[i].pos), btriangle[i].normal),0.0f);
	glm::vec4 s = glm::normalize(light - btriangle[i].pos) + glm::normalize(glm::vec4(10.0f,10.0f,10.0f,1.0f) - btriangle[i].pos);
 	float specular_intensity = max(glm::dot(glm::normalize(s),btriangle[i].pos), 0.0f);
	float R = btriangle[i].color.x * ( diffuse_intensity * diffusex + specular_intensity * specularx);
	float G = btriangle[i].color.y * ( diffuse_intensity * diffusey + specular_intensity * speculary);
	float B = btriangle[i].color.z * ( diffuse_intensity * diffusez + specular_intensity * specularz);
	btriangle[i].color = glm::vec4(R,G,B,1.0f);
		}

//calculate top
	for(int i = 0; i<360*12; i++){
	float diffuse_intensity = max(glm::dot(glm::normalize(light - ttriangle[i].pos), ttriangle[i].normal),0.0f);
	glm::vec4 s = glm::normalize(light - ttriangle[i].pos) + glm::normalize(glm::vec4(10.0f,10.0f,10.0f,1.0f) - ttriangle[i].pos);
 	float specular_intensity = max(glm::dot(glm::normalize(s),ttriangle[i].pos), 0.0f);
	float R = ttriangle[i].color.x * ( diffuse_intensity * diffusex + specular_intensity * specularx);
	float G = ttriangle[i].color.y * ( diffuse_intensity * diffusey + specular_intensity * speculary);
	float B = ttriangle[i].color.z * ( diffuse_intensity * diffusez + specular_intensity * specularz);
	ttriangle[i].color = glm::vec4(R,G,B,1.0f);
		}

//calculate first arm	
	for(int i = 0; i<36; i++){
	float diffuse_intensity = max(glm::dot(glm::normalize(light - ftriangle[i].pos), ftriangle[i].normal),0.0f);
	glm::vec4 s = glm::normalize(light - ftriangle[i].pos) + glm::normalize(glm::vec4(10.0f,10.0f,10.0f,1.0f) - ftriangle[i].pos);
 	float specular_intensity = max(glm::dot(glm::normalize(s),ftriangle[i].pos), 0.0f);
	float R = ftriangle[i].color.x * (diffuse_intensity * diffusex + specular_intensity * specularx);
	float G = ftriangle[i].color.y * (diffuse_intensity * diffusey + specular_intensity * speculary);
	float B = ftriangle[i].color.z * (diffuse_intensity * diffusez + specular_intensity * specularz);
	ftriangle[i].color = glm::vec4(R,G,B,1.0f);
		}

//calculate joint
	for(int i = 0; i<360*12; i++){
	float diffuse_intensity = max(glm::dot(glm::normalize(light - jtriangle[i].pos), jtriangle[i].normal),0.0f);
	glm::vec4 s = glm::normalize(light - jtriangle[i].pos) + glm::normalize(glm::vec4(10.0f,10.0f,10.0f,1.0f) - jtriangle[i].pos);
 	float specular_intensity = max(glm::dot(glm::normalize(s),jtriangle[i].pos), 0.0f);
	float R = jtriangle[i].color.x * ( diffuse_intensity * diffusex + specular_intensity * specularx);
	float G = jtriangle[i].color.y * ( diffuse_intensity * diffusey + specular_intensity * speculary);
	float B = jtriangle[i].color.z * ( diffuse_intensity * diffusez + specular_intensity * specularz);
	jtriangle[i].color = glm::vec4(R,G,B,1.0f);
		}

//calculate second arm
	for(int i = 0; i<360*12; i++){
	float diffuse_intensity = max(glm::dot(glm::normalize(light - striangle[i].pos), striangle[i].normal),0.0f);
	glm::vec4 s = glm::normalize(light - jtriangle[i].pos) + glm::normalize(glm::vec4(10.0f,10.0f,10.0f,1.0f) - striangle[i].pos);
 	float specular_intensity = max(glm::dot(glm::normalize(s),striangle[i].pos), 0.0f);
	float R = striangle[i].color.x * ( diffuse_intensity * diffusex + specular_intensity * specularx);
	float G = striangle[i].color.y * ( diffuse_intensity * diffusey + specular_intensity * speculary);
	float B = striangle[i].color.z * ( diffuse_intensity * diffusez + specular_intensity * specularz);
	striangle[i].color = glm::vec4(R,G,B,1.0f);
		}

//calculate pen
	for(int i = 0; i<360*12; i++){
	float diffuse_intensity = max(glm::dot(glm::normalize(light - ptriangle[i].pos), ptriangle[i].normal),0.0f);
	glm::vec4 s = glm::normalize(light - jtriangle[i].pos) + glm::normalize(glm::vec4(10.0f,10.0f,10.0f,1.0f) - ptriangle[i].pos);
 	float specular_intensity = max(glm::dot(glm::normalize(s),ptriangle[i].pos), 0.0f);
	float R = ptriangle[i].color.x * ( diffuse_intensity * diffusex + specular_intensity * specularx);
	float G = ptriangle[i].color.y * ( diffuse_intensity * diffusey + specular_intensity * speculary);
	float B = ptriangle[i].color.z * ( diffuse_intensity * diffusez + specular_intensity * specularz);
	ptriangle[i].color = glm::vec4(R,G,B,1.0f);
		}

//calculate box
	for(int i = 0; i<36; i++){
	float diffuse_intensity = max(glm::dot(glm::normalize(light - xtriangle[i].pos), xtriangle[i].normal),0.0f);
	glm::vec4 s = glm::normalize(light - xtriangle[i].pos) + glm::normalize(glm::vec4(10.0f,10.0f,10.0f,1.0f) - xtriangle[i].pos);
 	float specular_intensity = max(glm::dot(glm::normalize(s),xtriangle[i].pos), 0.0f);
	float R = xtriangle[i].color.x * ( diffuse_intensity * diffusex + specular_intensity * specularx);
	float G = xtriangle[i].color.y * ( diffuse_intensity * diffusey + specular_intensity * speculary);
	float B = xtriangle[i].color.z * ( diffuse_intensity * diffusez + specular_intensity * specularz);
	xtriangle[i].color = glm::vec4(R,G,B,1.0f);
		}

}


void initial_color(){
    for(int i=0; i<36; i++){
		ftriangle[i].color = glm::vec4(1,1,0.1,1);		
		xtriangle[i].color = glm::vec4(1,1,0.1,1);
		}
    for(int i=0; i<360*12; i++){
		btriangle[i].color = glm::vec4(0.5,0.5,0.5,1);		
		ttriangle[i].color = glm::vec4(0.1,0.1,1,1);
		jtriangle[i].color = glm::vec4(0.1,1,1,1);
		striangle[i].color = glm::vec4(1,0.1,1,1);
		ptriangle[i].color = glm::vec4(0.5,1,0.5,1);
		}
	}

void draw_grid(){
    
    GLuint vao;
    GLuint grid_vbo;
    
    GLuint position_location;
    GLuint color_location;
    GLuint normal_location;
    GLuint MVP_location;
    
    GLuint btriangle_vbo;
    GLuint ttriangle_vbo;
    GLuint ftriangle_vbo;
    GLuint jtriangle_vbo;
    GLuint striangle_vbo;
    GLuint ptriangle_vbo;
    GLuint xtriangle_vbo;
    GLuint htriangle_vbo;
    GLuint botriangle_vbo;
    GLuint atriangle1_vbo;
    GLuint atriangle2_vbo;
    GLuint ltriangle1_vbo;
    GLuint ltriangle2_vbo;

    float rotate_x; float rotate_y; float rotate_z;
    rotate_x=-cos(cx)*sin(cy)*10.0f;
    rotate_y=sin(cx)*10.0f;
    rotate_z=cos(cx)*cos(cy)*10.0f;
    

//add
//the camera is sitting at (10,10,10) and looking at (0,0,0), the up direction is (0,1,0)
    if(hint==0)
    {
    POS=glm::vec3(5.0f,5.0f,5.0f);
    }
    else if(hint==1)
    {
    POS=glm::vec3(rotate_x,rotate_y,rotate_z);
    }
    glm::vec3 AT=glm::vec3(0.0f,0.0f,0.0f);
    glm::vec3 UP=glm::vec3(0.0f,1.0f,0.0f);
    glm::mat4 view = glm::lookAt(POS,AT,UP);
//the view cube is in range [-5,5] on both horizontal and vertical direction, and [0.1, 100] on depth direction. 
    glm::mat4 projection=glm::ortho(-5.0f,5.0f,-5.0f,5.0f,0.1f,100.0f);
    glm::mat4 model =glm::mat4(1.0f);
    glm::mat4 MVP = projection*view*model;
       

    // specify the shaders we want to use
    glUseProgram(program);

    // get the input variable location in this shader
    position_location = glGetAttribLocation(program, "in_position");
    color_location = glGetAttribLocation(program, "in_color");
//add
    normal_location =glGetAttribLocation(program, "in_normal");
//add
    MVP_location = glGetUniformLocation(program, "MVP");    

    // create and bind a VAO
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // enable the input locations we wan to use
    glEnableVertexAttribArray(position_location);
    glEnableVertexAttribArray(color_location);
//    glEnableVertexAttribArray(normal_location);

//add   
    glUniformMatrix4fv(MVP_location, 1, GL_FALSE, &MVP[0][0]);

    // draw points

    // generate and bind a vertex buffer to hold the vertex data on GPU
    glGenBuffers(1, &grid_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, grid_vbo);

    // initialize the vertex buffer with the vertex data
    
    glBufferData(GL_ARRAY_BUFFER, grid_vertices.size() * sizeof(Vertex), &grid_vertices[0] , GL_STATIC_DRAW);

    // specify the input data format
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
//    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,normal));

    // draw points
    glPointSize(10);
    glDrawArrays(GL_LINES, 0, grid_vertices.size());


//add draw base
    glGenBuffers(1, &btriangle_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, btriangle_vbo);
    glBufferData(GL_ARRAY_BUFFER, 360*6*2*sizeof(Vertex), &btriangle[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
//    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));	
    glDrawArrays(GL_TRIANGLES, 0, 360*6*2);


//add draw top
    glGenBuffers(1, &ttriangle_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, ttriangle_vbo);
    glBufferData(GL_ARRAY_BUFFER, 360*6*2*sizeof(Vertex), &ttriangle[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
//    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));	
    glDrawArrays(GL_TRIANGLES, 0, 360*6*2);

//add draw firstarm
    glGenBuffers(1, &ftriangle_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, ftriangle_vbo);
    glBufferData(GL_ARRAY_BUFFER, 6*2*3*sizeof(Vertex), &ftriangle[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
//    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));	
    glDrawArrays(GL_TRIANGLES, 0, 6*2*3);

//add draw joint
    glGenBuffers(1, &jtriangle_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, jtriangle_vbo);
    glBufferData(GL_ARRAY_BUFFER, 360*6*2*sizeof(Vertex), &jtriangle[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
//    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));	
    glDrawArrays(GL_TRIANGLES, 0, 360*6*2);

//add draw secondarm
    glGenBuffers(1, &striangle_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, striangle_vbo);
    glBufferData(GL_ARRAY_BUFFER, 360*6*2*sizeof(Vertex), &striangle[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
//    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));	
    glDrawArrays(GL_TRIANGLES, 0, 360*6*2);

//add draw pen
    glGenBuffers(1, &ptriangle_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, ptriangle_vbo);
    glBufferData(GL_ARRAY_BUFFER, 360*6*2*sizeof(Vertex), &ptriangle[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
//    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));	
    glDrawArrays(GL_TRIANGLES, 0, 360*6*2);

//add draw box
    glGenBuffers(1, &xtriangle_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, xtriangle_vbo);
    glBufferData(GL_ARRAY_BUFFER, 6*2*3*sizeof(Vertex), &xtriangle[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
//    glVertexAttribPointer(normal_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));	
    glDrawArrays(GL_TRIANGLES, 0, 6*2*3);


//add draw head
    glGenBuffers(1, &htriangle_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, htriangle_vbo);
    glBufferData(GL_ARRAY_BUFFER, 6*2*3*sizeof(Vertex), &htriangle[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));	
    glDrawArrays(GL_TRIANGLES, 0, 6*2*3);

//add draw body
    glGenBuffers(1, &botriangle_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, botriangle_vbo);
    glBufferData(GL_ARRAY_BUFFER, 6*2*3*sizeof(Vertex), &botriangle[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
    glDrawArrays(GL_TRIANGLES, 0, 6*2*3);


//add draw leg1
    glGenBuffers(1, &ltriangle1_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, ltriangle1_vbo);
    glBufferData(GL_ARRAY_BUFFER, 6*2*3*sizeof(Vertex), &ltriangle1[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
    glDrawArrays(GL_TRIANGLES, 0, 6*2*3);

//add draw leg2
    glGenBuffers(1, &ltriangle2_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, ltriangle2_vbo);
    glBufferData(GL_ARRAY_BUFFER, 6*2*3*sizeof(Vertex), &ltriangle2[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
    glDrawArrays(GL_TRIANGLES, 0, 6*2*3);

//add draw arm1
    glGenBuffers(1, &atriangle1_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, atriangle1_vbo);
    glBufferData(GL_ARRAY_BUFFER, 6*2*3*sizeof(Vertex), &atriangle1[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
    glDrawArrays(GL_TRIANGLES, 0, 6*2*3);

//add draw arm2
    glGenBuffers(1, &atriangle2_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, atriangle2_vbo);
    glBufferData(GL_ARRAY_BUFFER, 6*2*3*sizeof(Vertex), &atriangle2[0], GL_STATIC_DRAW);
    glVertexAttribPointer(position_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), NULL);
    glVertexAttribPointer(color_location, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex,color));
    glDrawArrays(GL_TRIANGLES, 0, 6*2*3);


    // unbind VAO and VBO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // Delete VAO and VBO
    glDeleteBuffers(1, &grid_vbo);
    glDeleteBuffers(1, &btriangle_vbo);
    glDeleteBuffers(1, &ttriangle_vbo);
    glDeleteBuffers(1, &ftriangle_vbo);
    glDeleteBuffers(1, &jtriangle_vbo);
    glDeleteBuffers(1, &striangle_vbo);
    glDeleteBuffers(1, &ptriangle_vbo);
    glDeleteBuffers(1, &xtriangle_vbo);
    glDeleteBuffers(1, &htriangle_vbo);
    glDeleteBuffers(1, &botriangle_vbo);
    glDeleteBuffers(1, &atriangle1_vbo);
    glDeleteBuffers(1, &atriangle2_vbo);
    glDeleteBuffers(1, &ltriangle1_vbo);
    glDeleteBuffers(1, &ltriangle2_vbo);
    glDeleteVertexArrays(1, &vao);
  
}


void display(){
    // Clear Viewport
    glClearColor(0.0f,0.0f,0.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if(light==1){
    lighting();
    }
    else if(light==0)
    {initial_color();
    }

    draw_grid();

    glFlush();
    glutSwapBuffers();

}

void reshape(int width, int height){
    // Clip the view port to match our ratio
    glViewport(0, 0, width, height);
    glutPostRedisplay();
}

void graphics_init(){

    // init vertex shader
    read_shader_source_code("vs.glsl", buf, 1024);
    cout << buf << endl;
    int vertex_shader_source_length = strlen(buf);
    const char *vertex_shader_sources[] = { buf };

    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_sources, &vertex_shader_source_length);
    glCompileShader(vertex_shader);
    
    // init fragment shader 
    read_shader_source_code("fs.glsl", buf, 1024); 
    int fragment_shader_source_length = strlen(buf);
    const char *fragment_shader_sources[] = { buf };
    
    fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, fragment_shader_sources, &fragment_shader_source_length);
    glCompileShader(fragment_shader);

    // init program
    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);

    glLinkProgram(program);
  
    // enable depth test
    glEnable(GL_DEPTH_TEST);
    
    glEnable(GL_LIGHTING);
}

void print_info()
{
    fprintf(
        stdout,
        "INFO: OpenGL Version: %s\n",
        glGetString(GL_VERSION)
        );
    fprintf(stdout, "INFO: Using GLEW %s\n", glewGetString(GLEW_VERSION));

    if(glewIsSupported("GL_ARB_debug_output"))
    {
        fprintf(stdout, "INFO: Support ARB_debug_output\n");
    }
    else
    {
        fprintf(stdout, "INFO: Not support ARB_debug_output\n");
    }
}

int main(int argc,char * argv[]){

    // Setup GLUT
    glutInit(&argc,argv);

    glutInitContextVersion (3, 0);
    glutInitContextFlags (GLUT_DEBUG);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB     | GLUT_DEPTH | GLUT_STENCIL);
    glutInitWindowSize (window_width,window_height);
    glutCreateWindow("First step - OpenGL 3");
  
    // Setup GLEW
    glewExperimental = true;
    GLenum err = glewInit();
    if(GLEW_OK != err)
    {
        /* Problem: glewInit failed, something is seriously wrong. */
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
    }

    print_info();

    if(glewIsSupported("GL_ARB_debug_output"))
        VSDebugLib::init();
    // Initialize OpenGL
    graphics_init ();

    initialize_points();
    initialize_base();
    initialize_top();
    initialize_firstarm();
    initialize_joint();
    initialize_secondarm();
    initialize_pen();
    initialize_box();
    initialize_head();
    initialize_body();
    initialize_arm1();
    initialize_arm2();
    initialize_leg1();
    initialize_leg2();
    // set up callbacks
    glutReshapeFunc(&reshape);
    glutDisplayFunc(&display);
    glutKeyboardFunc(&keyboard);
    glutMouseFunc(&mouse);
    glutMotionFunc(&motion);
    glutSpecialFunc(&special_key);

    // main loop
    glutMainLoop();
    return 0;
}
