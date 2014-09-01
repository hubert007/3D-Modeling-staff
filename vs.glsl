#version 130

precision mediump float;

in vec4 in_position;
in vec4 in_color;
in vec4 in_normal;

out vec4 vColor;

uniform mat4 MVP;

void main()
{
 
    gl_Position = MVP * in_position;
    vColor = in_color;
}
