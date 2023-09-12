
#include <GL/glut.h>
#include <bits/stdc++.h>
#include "bitmap_image.hpp"
#include "1805065_classes.cpp"

#define pi (2 * acos(0.0))

using namespace std;

double angle = 0;
int windowWidth = 500;
int windowHeight = 500;
point pos; // position of the eye 0 -100 20
point l;   // look/forward direction 0 1 0
point r;   // right direction
point u;   // up direction
point checkerBoardPos;
double fovY, aspectRatio;
double nearDist, farDist;
vector<struct point> pointBuffer;
vector<Object *> objects;
vector<Light *> lights;
vector<SpotLight *> spotLights;
int levelOfRecursion, numberOfPixels;
bitmap_image image;
bitmap_image textureImageW, textureImageB;
int loadTexture = 0;

void loadData()
{
    textureImageW = bitmap_image("./texture_w.bmp");
    textureImageB = bitmap_image("./texture_b.bmp");
    ifstream sceneFile("./description.txt");
    sceneFile >> nearDist >> farDist;
    sceneFile >> fovY;
    sceneFile >> aspectRatio;
    sceneFile >> levelOfRecursion;
    sceneFile >> numberOfPixels;

    double color[3], coEfficients[4];
    int shine;
    double specularCoEfficient;
    double t1, t2;
    point p1, p2;
    // First object is checker board
    // Get width of the checker board
    sceneFile >> t1;
    sceneFile >> coEfficients[0] >> coEfficients[1] >> coEfficients[2];
    // cout << "CheckerBoard: " << t1 << coEfficients[0] << coEfficients[1] << coEfficients[2] << endl;
    CheckerBoard *checkerBoard = new CheckerBoard(t1, coEfficients);
    objects.push_back(checkerBoard);
    sceneFile >> t1;
    // cout << *checkerBoard << endl;
    int numberOfObjects = (int)t1;

    for (int i = 0; i < numberOfObjects; i++)
    {
        string objectType;
        sceneFile >> objectType;

        if (objectType == "sphere")
        {
            sceneFile >> p1.x >> p1.y >> p1.z;
            sceneFile >> t1;
            sceneFile >> color[0] >> color[1] >> color[2];
            sceneFile >> coEfficients[0] >> coEfficients[1] >> specularCoEfficient >> coEfficients[2];
            sceneFile >> shine;
            coEfficients[3] = coEfficients[2];
            coEfficients[2] = specularCoEfficient;
            Sphere *sphere = new Sphere(p1, t1, color, coEfficients, shine);
            objects.push_back(sphere);
        }
        else if (objectType == "pyramid")
        {
            sceneFile >> p1.x >> p1.y >> p1.z;
            sceneFile >> t1 >> t2;
            sceneFile >> color[0] >> color[1] >> color[2];
            sceneFile >> coEfficients[0] >> coEfficients[1] >> specularCoEfficient >> coEfficients[2];
            sceneFile >> shine;
            coEfficients[3] = coEfficients[2];
            coEfficients[2] = specularCoEfficient;
            Pyramid *pyramid = new Pyramid(p1, t1, t2, color, coEfficients, shine);
            objects.push_back(pyramid);
        }
        else if (objectType == "cube")
        {
            sceneFile >> p1.x >> p1.y >> p1.z;
            sceneFile >> t1;
            sceneFile >> color[0] >> color[1] >> color[2];
            sceneFile >> coEfficients[0] >> coEfficients[1] >> specularCoEfficient >> coEfficients[2];
            sceneFile >> shine;
            coEfficients[3] = coEfficients[2];
            coEfficients[2] = specularCoEfficient;
            Cube *cube = new Cube(p1, t1, color, coEfficients, shine);
            objects.push_back(cube);
        }
    }

    sceneFile >> t1;
    int numberOfLights = (int)t1;
    for (int i = 0; i < numberOfLights; i++)
    {
        sceneFile >> p1.x >> p1.y >> p1.z;
        sceneFile >> t2;
        Light *light = new Light(p1, t2);
        lights.push_back(light);
    }

    sceneFile >> t1;
    numberOfLights = (int)t1;
    for (int i = 0; i < numberOfLights; i++)
    {
        sceneFile >> p1.x >> p1.y >> p1.z;
        sceneFile >> t2;
        sceneFile >> p2.x >> p2.y >> p2.z;
        sceneFile >> t1;
        SpotLight *spotLight = new SpotLight(p1, t2, p2, t1);
        spotLights.push_back(spotLight);
    }
}

pair<double, double> getLambertAndPhongCoefficients(point intersect,int nearest)
{
    double lambert = 0, phong = 0;
    // cout << "lights.size(): " << lights.size() << endl;
    for (int i = 0; i < lights.size(); i++)
    {
        point lightDirection = intersect - lights[i]->pos ;
        lightDirection.normalize();

        Ray ray(lights[i]->pos, lightDirection);
        double t_curr = objects[nearest]->intersect(ray);
        double t = 0;
        bool illuminate = true;
        for(int j = 0; j<objects.size(); j++){
            if(j == nearest) continue;
            t = objects[j]->intersect(ray);
            if(t > 0 && t < t_curr) {
                illuminate = false;
                break;
            }
        }
        if(!illuminate) continue;

        point normal = objects[nearest]->getNormal(intersect);

        if (ray.dir * normal > 0)
            normal = normal * -1;
        double distance = intersect.distance(lights[i]->pos);

        double scalingfactor = exp(-(distance * distance * lights[i]->falloff));
        point toSource = lights[i]->pos - intersect;
        toSource.normalize();
        double lambertCoefficient = max(0.0, normal * toSource);

        lambert += lambertCoefficient * scalingfactor;

        // point reflectionDirection = toSource - normal * (2 * (toSource * normal));
        point reflectionDirection = normal * (2 * (toSource * normal)) - toSource;
        reflectionDirection.normalize();
        double phongCoefficient = max(0.0, reflectionDirection * toSource);
        phong += pow(phongCoefficient, objects[nearest]->shine) * scalingfactor;
    }

    for (int i = 0; i < spotLights.size(); i++)
    {
        point lightDirection = intersect - spotLights[i]->pos;
        lightDirection.normalize();

        spotLights[i]->dir.normalize();
        double angle = acos(lightDirection * spotLights[i]->dir) * 180 / pi;

        if(angle > spotLights[i]->cutoff) continue;

        Ray ray(spotLights[i]->pos, lightDirection);
        double t_curr = objects[nearest]->intersect(ray);
        // cout << "t_curr: " << t_curr << endl;
        double t = 0;
        bool illuminate = true;
        for(int j = 0; j<objects.size(); j++){
            if(j == nearest) continue;
            t = objects[j]->intersect(ray);
            if(t > 0 && t < t_curr) {
                illuminate = false;
                break;
            }
        }
        if(!illuminate) continue;

        point normal = objects[nearest]->getNormal(intersect);
        if (ray.dir * normal > 0)
            normal = normal * -1;

        double distance = intersect.distance(spotLights[i]->pos);

        double scalingfactor = exp(-(distance * distance * spotLights[i]->falloff));
        point toSource = spotLights[i]->pos - intersect;
        toSource.normalize();
        double lambertCoefficient = max(0.0, normal * toSource);
        lambert += lambertCoefficient * scalingfactor;

        // point reflectionDirection = toSource - normal * (2 * (toSource * normal));
        point reflectionDirection = normal * (2 * (toSource * normal)) - toSource;
        reflectionDirection.normalize();
        double phongCoefficient = max(0.0, reflectionDirection * toSource);
        phong += pow(phongCoefficient, objects[nearest]->shine) * scalingfactor;
        
    }
    return make_pair(lambert, phong);
}

point rayTrace(Ray ray,int level){
    if(level <=0 ) return point(0,0,0);
    double t, t_min = 10000;
    int nearest = -1;
    for (int k = 0; k < objects.size(); k++)
    {
        t = objects[k]->intersect(ray);
        if (t > 0 && t < t_min)
        {
            t_min = t;
            nearest = k;
        }
    }
    if (nearest == -1)
    {
        return point(0,0,0);
    }

    double phong = 0 , lambert = 0;
    point intersectionPoint = ray.getPoint(t_min);
    point normal = objects[nearest]->getNormal(intersectionPoint);
    //Get lambert and phong coefficients
    pair<double, double> lambertAndPhongCoefficients = getLambertAndPhongCoefficients(intersectionPoint, nearest);
    lambert = lambertAndPhongCoefficients.first;
    phong = lambertAndPhongCoefficients.second;

    point color = objects[nearest]->getColor(ray.getPoint(t_min));
    double coefficients= objects[nearest]->coEfficients[0] + objects[nearest]->coEfficients[1] * lambert + objects[nearest]->coEfficients[2] * phong;
    // image.set_pixel(i, j, color.x * 255, color.y * 255, color.z * 255);
    point reflectionDirection = ray.dir - normal * (2 * (ray.dir * normal));
    reflectionDirection.normalize();
    Ray reflectionRay(intersectionPoint, reflectionDirection);
    point reflectionColor = rayTrace(reflectionRay, level - 1);
    return color * coefficients + reflectionColor * objects[nearest]->coEfficients[3];
    
}

void Capture()
{
    cout << "Capturing" << endl;
    cout << "fovY: " << fovY << endl;
    // initPointBuffer();

    image = bitmap_image(numberOfPixels, numberOfPixels);

    double height = 2 * nearDist * tan(fovY / 2 * pi / 180);
    double width = height * aspectRatio;
    double dy = height / numberOfPixels;
    double dx = width / numberOfPixels;

    point topLeft = pos + l * nearDist + u * (height / 2) - r * (width / 2);
    topLeft = topLeft + r * dx / 2 - u * dy / 2;

    for(int  i=0;i< numberOfPixels;i++){
        for(int j=0;j<numberOfPixels;j++){
            
            point currentPixel = topLeft + r * dx * i - u * dy * j;
            Ray ray(pos, currentPixel - pos);

            point color = rayTrace(ray, levelOfRecursion);
            image.set_pixel(i, j, color.x * 255, color.y * 255, color.z * 255);
        }
    }

    for (int i = 0; i < lights.size(); i++)
        cout << *lights[i] << endl;

    for (int i = 0; i < spotLights.size(); i++)
        cout << *spotLights[i] << endl;
    image.save_image("./1805065_output.bmp");
    cout << "Captured" << endl;
}

void keyboardListener(unsigned char key, int xx, int yy)
{
    double rate = 15 * pi / 180;
    point tempL = l;
    point tempR = r;
    point tempU = u;
    tempL.normalize();
    tempR.normalize();
    tempU.normalize();
    switch (key)
    {
    case '0':
        Capture();
        break;

    case ' ':
        loadTexture = 1 - loadTexture;
        cout << "loadTexture: " << loadTexture << endl;
        break;

    case '1':
        r = r * cos(rate) + l * sin(rate);
        l = l * cos(rate) - tempR * sin(rate);
        // r.x = r.x * cos(rate) + l.x * sin(rate);
        // r.y = r.y * cos(rate) + l.y * sin(rate);
        // r.z = r.z * cos(rate) + l.z * sin(rate);

        // l.x = l.x * cos(rate) - r.x * sin(rate);
        // l.y = l.y * cos(rate) - r.y * sin(rate);
        // l.z = l.z * cos(rate) - r.z * sin(rate);
        break;

    case '2':
        r = r * cos(-rate) + l * sin(-rate);
        l = l * cos(-rate) - tempR * sin(-rate);
        // r.x = r.x * cos(-rate) + l.x * sin(-rate);
        // r.y = r.y * cos(-rate) + l.y * sin(-rate);
        // r.z = r.z * cos(-rate) + l.z * sin(-rate);

        // l.x = l.x * cos(-rate) - r.x * sin(-rate);
        // l.y = l.y * cos(-rate) - r.y * sin(-rate);
        // l.z = l.z * cos(-rate) - r.z * sin(-rate);
        break;

    case '3':
        l = l * cos(rate) + u * sin(rate);
        u = u * cos(rate) - tempL * sin(rate);
        // l.x = l.x * cos(rate) + u.x * sin(rate);
        // l.y = l.y * cos(rate) + u.y * sin(rate);
        // l.z = l.z * cos(rate) + u.z * sin(rate);

        // u.x = u.x * cos(rate) - l.x * sin(rate);
        // u.y = u.y * cos(rate) - l.y * sin(rate);
        // u.z = u.z * cos(rate) - l.z * sin(rate);
        break;

    case '4':
        l = l * cos(-rate) + u * sin(-rate);
        u = u * cos(-rate) - tempL * sin(-rate);
        // l.x = l.x * cos(-rate) + u.x * sin(-rate);
        // l.y = l.y * cos(-rate) + u.y * sin(-rate);
        // l.z = l.z * cos(-rate) + u.z * sin(-rate);

        // u.x = u.x * cos(-rate) - l.x * sin(-rate);
        // u.y = u.y * cos(-rate) - l.y * sin(-rate);
        // u.z = u.z * cos(-rate) - l.z * sin(-rate);
        break;

    case '5':
        u = u * cos(rate) + r * sin(rate);
        r = r * cos(rate) - tempU * sin(rate);
        // u.x = u.x * cos(rate) + r.x * sin(rate);
        // u.y = u.y * cos(rate) + r.y * sin(rate);
        // u.z = u.z * cos(rate) + r.z * sin(rate);

        // r.x = r.x * cos(rate) - u.x * sin(rate);
        // r.y = r.y * cos(rate) - u.y * sin(rate);
        // r.z = r.z * cos(rate) - u.z * sin(rate);
        break;

    case '6':
        u = u * cos(-rate) + r * sin(-rate);
        r = r * cos(-rate) - tempU * sin(-rate);
        // u.x = u.x * cos(-rate) + r.x * sin(-rate);
        // u.y = u.y * cos(-rate) + r.y * sin(-rate);
        // u.z = u.z * cos(-rate) + r.z * sin(-rate);

        // r.x = r.x * cos(-rate) - u.x * sin(-rate);
        // r.y = r.y * cos(-rate) - u.y * sin(-rate);
        // r.z = r.z * cos(-rate) - u.z * sin(-rate);
        break;

    default:
        break;
    }
    l.normalize();
    r.normalize();
    u.normalize();
    glutPostRedisplay();
}

void specialKeyListener(int key, int x, int y)
{
    double step = 10;
    switch (key)
    {
    case GLUT_KEY_UP: // down arrow key
        pos = pos + l * step;
        break;
    case GLUT_KEY_DOWN: // up arrow key
        pos = pos - l * step;
        break;

    case GLUT_KEY_RIGHT:
        pos = pos + r * step;
        break;
    case GLUT_KEY_LEFT:
        pos = pos - r * step;
        break;

    case GLUT_KEY_PAGE_UP:
        pos = pos + u * step;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos = pos - u * step;
        break;

    default:
        break;
    }
    glutPostRedisplay();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0); // color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    r.x = l.y * u.z - l.z * u.y;
    r.y = l.z * u.x - l.x * u.z;
    r.z = l.x * u.y - l.y * u.x;

    gluLookAt(pos.x, pos.y, pos.z,
              pos.x + l.x, pos.y + l.y, pos.z + l.z,
              u.x, u.y, u.z);

    glMatrixMode(GL_MODELVIEW);

    for(int i = 0; i < objects.size(); i++) {
    	objects[i]->draw();
    }


    // //Test ray
    // point start = {0, 0, 0};
    // point dir = {1, 1, 1};
    // Ray ray = {start, dir};
    // cout << ray << endl;

    // //Draw ray
    // glBegin(GL_LINES);
    // 	glColor3f(1, 0, 0);
    // 	glVertex3f(ray.start.x, ray.start.y, ray.start.z);
    // 	glVertex3f(ray.start.x + ray.dir.x, ray.start.y + ray.dir.y, ray.start.z + ray.dir.z);
    // glEnd();

    // drawAxes(50);
    // drawgrid=1;
    // int width = 50;
    // drawGrid(width);

    // point center = {20,20,20};
    // int radius = 20;
    // point cubeBottomLowerLeft = {-100, -100, 10};
    // double side = 40;
    // glPushMatrix();
    // {
    // 	glTranslatef(cubeBottomLowerLeft.x, cubeBottomLowerLeft.y, cubeBottomLowerLeft.z);
    // 	glScalef(side, side, side);
    // 	drawCube(cubeBottomLowerLeft, side, 0.0, 0.5, 1.0);
    // }
    // glPopMatrix();
    // drawSphere(center, radius, 20, 20, 0.25, 0.3, 1.0);

    // center = {-20.0, -20.0, 20.0};
    // radius = 15;
    // drawSphere(center, radius, 20, 20, 1.0, 0.0, 1.0);
    // glPushMatrix();
    // {
    // 	point lowerLeft = {-40.0, 0.0, 5.0};
    // 	double width = 30.0;
    // 	double height = 40.0;
    // 	glTranslatef(lowerLeft.x, lowerLeft.y, lowerLeft.z);
    // 	glScalef(width, width, height);
    // 	drawPyramid(1.0, 0.0, 0.0);
    // }
    // glPopMatrix();
    glutSwapBuffers();
}

void animate()
{
    angle += 0.05;
    // codes for any changes in Models, Camera
    glutPostRedisplay();
}

void initGL()
{
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); 
    glEnable(GL_DEPTH_TEST); // Enable depth testing for z-culling
}

void reshape(GLsizei width, GLsizei height)
{ // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0)
        height = 1; // To prevent divide by 0
    // GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity();
    gluPerspective(fovY, aspectRatio, nearDist, farDist);
}

int main(int argc, char **argv)
{
    loadData();

    // pos = {0, -100, 20};
    // l = {0, 1, 0};
    // r = {1, 0, 0};
    // u = {0, 0, 1};

    pos = point(0, -100, 20);
    checkerBoardPos = point(0, -100, 20);
    l = point(0, 1, 0);
    r = point(1, 0, 0);
    u = point(0, 0, 1);

    glutInit(&argc, argv); // Initialize GLUT
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(50, 50);                           // Position the window's initial top-left corner
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color
    glutCreateWindow("RayTracing - 1805065");
    glutDisplayFunc(display); // Register display callback handler for window re-paint
    glutReshapeFunc(reshape); // Register callback handler for window re-shape

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    // cout<< "Main loop ended " << numberOfPixels;
    initGL();       // Our own OpenGL initialization
    glutMainLoop(); // Enter the event-processing loop

    objects.clear();
    objects.shrink_to_fit();

    return 0;
}