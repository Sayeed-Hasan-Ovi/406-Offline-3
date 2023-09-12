#include <GL/glut.h>
#include <bits/stdc++.h>
#include "bitmap_image.hpp"
#include "1805065_classes.cpp"
#include "termcolor.hpp"

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
vector<Object *> objects;
vector<Light *> lights;
vector<SpotLight *> spotLights;
int levelOfRecursion, numberOfPixels;
bitmap_image image;
bitmap_image textureImageW, textureImageB;
int loadTexture = 0;

void colorify(string text, string color)
{
    if (color == "red")
    {
        cout << termcolor::bold << termcolor::red;
    }
    else if (color == "green")
    {
        cout << termcolor::bold << termcolor::green;
    }
    else if (color == "blue")
    {
        cout << termcolor::bold << termcolor::blue;
    }
    else if (color == "yellow")
    {
        cout << termcolor::bold << termcolor::yellow;
    }
    else if (color == "cyan")
    {
        cout << termcolor::bold << termcolor::cyan;
    }
    else if (color == "magenta")
    {
        cout << termcolor::bold << termcolor::magenta;
    }
    else if (color == "white")
    {
        cout << termcolor::bold << termcolor::white;
    }
    else if (color == "grey")
    {
        cout << termcolor::bold << termcolor::grey;
    }
    else if (color == "reset")
    {
        cout << termcolor::reset;
    }
    cout << text << endl;
    cout << termcolor::reset;
}

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

    sceneFile >> t1; // Size of checkerboard
    sceneFile >> coEfficients[0] >> coEfficients[1] >> coEfficients[2];

    CheckerBoard *checkerBoard = new CheckerBoard(t1, coEfficients);
    objects.push_back(checkerBoard);
    sceneFile >> t1;

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
    sceneFile.close();
    colorify("Data Loaded Successfully!", "green");
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
    // CAPTURING with fancy text
    cout << endl
         << endl
         << endl
         << endl
         << endl;
    string text = "                    Capturing...        ";
    colorify(text, "red");
    float progress{0.0f};
    size_t bar_width{60};
    string fill{"#"}, remainder{" "}, status_text{""};

    auto set_progress = [&](float value)
    {
        progress = value;
    };

    auto set_bar_width = [&](size_t width)
    {
        bar_width = width;
    };

    auto fill_bar_progress_with = [&](const string &chars)
    {
        fill = chars;
    };

    auto fill_bar_remainder_with = [&](const string &chars)
    {
        remainder = chars;
    };

    auto set_status_text = [&](const string &status)
    {
        status_text = status;
    };

    auto update = [&](float value, ostream &os = cout)
    {
        set_progress(value);
        os << '\r' << termcolor::bold << termcolor::green << '[';
        size_t pos = size_t(bar_width * progress / 100.0f);
        for (size_t i = 0; i < bar_width; ++i)
        {
            if (i < pos)
                os << fill;
            else if (i == pos)
                os << termcolor::reset << termcolor::bold << termcolor::red << ">";
            else
                os << remainder;
        }
        os << ']' << termcolor::reset << ' ' << setprecision(4) << progress << "% " << status_text;
        os.flush();
    };

    image = bitmap_image(numberOfPixels, numberOfPixels);

    double height = 2 * nearDist * tan(fovY / 2 * pi / 180);
    double width = height * aspectRatio;
    double dy = height / numberOfPixels;
    double dx = width / numberOfPixels;

    point topLeft = pos + l * nearDist + u * (height / 2) - r * (width / 2);
    topLeft = topLeft + r * dx / 2 - u * dy / 2;

    set_status_text("Loading...");
    int total = numberOfPixels * numberOfPixels;
    int count = 0;

    for (int i = 0; i < numberOfPixels; i++)
    {
        for (int j = 0; j < numberOfPixels; j++)
        {
            if (i == 0 && j == 0)
            {
                count = 0;
            }
            else
            {
                count = (i * numberOfPixels + j) * 100.0 / total;
            }
            point currentPixel = topLeft + r * dx * i - u * dy * j;
            Ray ray(pos, currentPixel - pos);

            point color = rayTrace(ray, levelOfRecursion);
            image.set_pixel(i, j, color.x * 255, color.y * 255, color.z * 255);
            update(count);
        }
    }
    cout << endl;

    image.save_image("./1805065_output.bmp");

    string text2 = "                    Captured!              ";

    cout << endl
         << endl
         << endl
         << endl
         << endl;
    colorify(text2, "green");

    cout << endl
         << endl
         << endl
         << endl
         << endl;
}

void keyboardListener(unsigned char key, int xx, int yy)
{
    double rate = 5 * pi / 180;
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
        break;

    case '2':
        r = r * cos(-rate) + l * sin(-rate);
        l = l * cos(-rate) - tempR * sin(-rate);
        break;

    case '3':
        l = l * cos(rate) + u * sin(rate);
        u = u * cos(rate) - tempL * sin(rate);
        break;

    case '4':
        l = l * cos(-rate) + u * sin(-rate);
        u = u * cos(-rate) - tempL * sin(-rate);
        break;

    case '5':
        u = u * cos(rate) + r * sin(rate);
        r = r * cos(rate) - tempU * sin(rate);
        break;

    case '6':
        u = u * cos(-rate) + r * sin(-rate);
        r = r * cos(-rate) - tempU * sin(-rate);
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
    double step = 5;
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

    for (int i = 0; i < objects.size(); i++)
    {
        objects[i]->draw();
    }
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

    pos = point(0, -100, 60);
    checkerBoardPos = point(0, -100, 60);
    l = point(0, 1, 0);
    r = point(1, 0, 0);
    u = point(0, 0, 1);

    glutInit(&argc, argv);
    glutInitWindowSize(windowWidth, windowHeight);
    // glutInitWindowPosition(50, 50);                           // Position the window's initial top-left corner
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color
    glutCreateWindow("RayTracing - 1805065");
    glutDisplayFunc(display);
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