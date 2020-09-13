#include <vector>
#include <math.h>

//Trigonometric constants
const float sqrt3 = sqrt(3.0);
const float sin120 = sqrt3/2.0;   
const float cos120 = -0.5;        
const float tan60 = sqrt3;
const float sin30 = 0.5;
const float tan30 = 1/sqrt3;

//Robot dimensions
float effector_side;
float base_side;

float length_a;
float length_b;
float length_c;
float length_d;

float current_position_x;
float current_position_y;
float current_position_z;

float next_position_x;
float next_position_y;
float next_position_z;

void loop() {
    //Obliczenie kąta o jaki trzeba obrócić serwo żeby przesunąć efektor we wskazane miejsce
    float gamma_A = inverseKinematics(next_position_x, next_position_y, next_position_z, 'A');
    float gamma_B = inverseKinematics(next_position_x, next_position_y, next_position_z, 'B');
    float gamma_C = inverseKinematics(next_position_x, next_position_y, next_position_z, 'C');       
}

/* Wyjaśnienie obliczeń zadania odwrotnego kinematyki:

*/

float inverseKinematics(float posX, float posY, float posZ, char servo)
{
    float x = 0.0, y = 0.0, z = 0.0;
    float M_PI120 = 120.0 * (M_PI / 180.0);
    float M_PI240 = 240.0 * (M_PI / 180.0);

    if (servo == 'A'){
        x = posX;
        y = posY;
        z = posZ;
    }

    if (servo == 'B'){
        x =  (cos(M_PI120)*(posX)) + (sin(M_PI120)*(posY));
        y = -(sin(M_PI120)*(posX)) + (cos(M_PI120)*(posY));
        z = posZ;
    }

    if (servo == 'C'){
        x =  (cos(M_PI240)*(posX)) + (sin(M_PI240)*(posY));
        y = -(sin(M_PI240)*(posX)) + (cos(M_PI240)*(posY));
        z = posZ;
    }

    float length1 = (length_a - length_d - y);
    float alpha = (360.0 / (2.0 * M_PI)) * (atan2(z, length1));
    float length2 = sqrt(pow(length1, 2.0) + pow(z, 2.0));
    float lenght3 = sqrt(pow(length_c, 2.0) - pow(x, 2.0));
    float beta = (360.0 / (2.0 * M_PI)) * (acos((pow(lenght3, 2) - pow(length2, 2.0) - pow(length_b, 2.0)) / (-2.0 * length2 * length_b)));
    float gamma = 180.0 - alpha - beta;

    return gamma;
}

/* Wyjaśnienie obliczeń zadania prostego kinematyki:
Now the three joint angles theta1, theta2 and theta3 are given, and we need to find the coordinates (x0, y0, z0) of end effector point E0.
As we know angles theta, we can easily find coordinates of points J1, J2 and J3 (see fig. below). Joints J1E1, J2E2 and J3E3 can freely rotate around points J1, J2 and J3 respectively, forming three spheres with radius re.
Now let's do the following: move the centers of the spheres from points J1, J2 and J3 to the points J'1, J'2 and J'3 using transition vectors E1E0, E2E0 and E3E0 respectively. After this transition all three spheres will intersect in one point: E0, as it is shown in fig. below:
So, to find coordinates (x0, y0, z0) of point E0, we need to solve set of three equations like (x-xj)^2+(y-yj)^2+(z-zj)^2 = re^2, where coordinates of sphere centers (xj, yj, zj) and radius re are known.
First, let's find coordinates of points J'1, J'2, J'3:
In the following equations I'll designate coordinates of points J1, J2, J3 as (x1, y1, z1), (x2, y2, z2) and (x3, y3, z3). Please note that x0=0. Here are equations of three spheres:
Finally, we need to solve this quadric equation and find z0 (we should choose the smallest negative equation root), and then calculate x0 and y0 from eq. (7) and (8).
*/

int forwardKinematics(float theta1, float theta2, float theta3, float &x0, float &y0, float &z0) {
     float t = (base_side-effector_side)*tan30/2;
     float dtr = M_PI/(float)180.0;
 
     theta1 *= dtr;
     theta2 *= dtr;
     theta3 *= dtr;
 
     float y1 = -(t + length_b*cos(theta1));
     float z1 = -length_b*sin(theta1);
 
     float y2 = (t + length_b*cos(theta2))*sin30;
     float x2 = y2*tan60;
     float z2 = -length_b*sin(theta2);
 
     float y3 = (t + length_b*cos(theta3))*sin30;
     float x3 = -y3*tan60;
     float z3 = -length_b*sin(theta3);
 
     float dnm = (y2-y1)*x3-(y3-y1)*x2;
 
     float w1 = y1*y1 + z1*z1;
     float w2 = x2*x2 + y2*y2 + z2*z2;
     float w3 = x3*x3 + y3*y3 + z3*z3;
     
     // x = (a1*z + b1)/dnm
     float a1 = (z2-z1)*(y3-y1)-(z3-z1)*(y2-y1);
     float b1 = -((w2-w1)*(y3-y1)-(w3-w1)*(y2-y1))/2.0;
 
     // y = (a2*z + b2)/dnm;
     float a2 = -(z2-z1)*x3+(z3-z1)*x2;
     float b2 = ((w2-w1)*x3 - (w3-w1)*x2)/2.0;
 
     // a*z^2 + b*z + c = 0
     float a = a1*a1 + a2*a2 + dnm*dnm;
     float b = 2*(a1*b1 + a2*(b2-y1*dnm) - z1*dnm*dnm);
     float c = (b2-y1*dnm)*(b2-y1*dnm) + b1*b1 + dnm*dnm*(z1*z1 - length_c*length_c);
  
     // discriminant
     float d = b*b - (float)4.0*a*c;
     if (d < 0) return -1; // non-existing point
 
     z0 = -(float)0.5*(b+sqrt(d))/a;
     x0 = (a1*z0 + b1)/dnm;
     y0 = (a2*z0 + b2)/dnm;
     return 0;
 }