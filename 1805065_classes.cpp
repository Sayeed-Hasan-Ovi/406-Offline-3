
#include <GL/glut.h>
#include <bits/stdc++.h>
#include "bitmap_image.hpp"
#define pi (2 * acos(0.0))

using namespace std;
extern bitmap_image textureImageW;
extern bitmap_image textureImageB;
extern int loadTexture;

class point
{
public:
    double x, y, z;

    point() {
        x = 0;
        y = 0;
        z = 0;
    }

    point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    point operator=(const point &p)
    {
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }

    point operator+(const point &p) const
    {
        point temp;
        temp.x = x + p.x;
        temp.y = y + p.y;
        temp.z = z + p.z;
        return temp;
    }

    point operator-(const point &p) const
    {
        point temp;
        temp.x = x - p.x;
        temp.y = y - p.y;
        temp.z = z - p.z;
        return temp;
    }

    point operator*(double d)
    {
        point temp;
        temp.x = x * d;
        temp.y = y * d;
        temp.z = z * d;
        return temp;
    }

    point operator/(const double &d) const
    {
        point temp;
        temp.x = x / d;
        temp.y = y / d;
        temp.z = z / d;
        return temp;
    }

    // cross product
    point operator^(point p)
    {
        point temp;
        temp.x = y * p.z - z * p.y;
        temp.y = z * p.x - x * p.z;
        temp.z = x * p.y - y * p.x;
        return temp;
    }

    // dot product
    double operator*(point p)
    {
        return x * p.x + y * p.y + z * p.z;
    }

    void normalize()
    {
        double length = sqrt(x * x + y * y + z * z);
        if (length <= 1e-6)
            return;
        x /= length;
        y /= length;
        z /= length;
    }

    double distance(point p)
    {
        return sqrt((x - p.x) * (x - p.x) +
                    (y - p.y) * (y - p.y) +
                    (z - p.z) * (z - p.z));
    }
};

class Ray{
    public:
        point start;
        point dir;
        Ray(point start, point dir){
            this->start = start;
            this->dir = dir;
            this->dir.normalize();
        }

        point getPoint(double t){
            return start + dir * t;
        }
};

class Object{
    public:
        double coEfficients[4];
        point color;
        double shine;
        Object()
        {
        coEfficients[0] = 0; //ambient
        coEfficients[1] = 0; //diffuse
        coEfficients[2] = 0; //specular
        coEfficients[3] = 0; //reflection
        }
        Object(double coEfficients[4])
        {
        this->coEfficients[0] = coEfficients[0];
        this->coEfficients[1] = coEfficients[1];
        this->coEfficients[2] = coEfficients[2];
        this->coEfficients[3] = coEfficients[3];
        }
        virtual void draw() = 0;
        // virtual double intersect(Ray ray, color &finalColor, int levelOfRecursion) = 0;
        virtual point getColor(point p) = 0; 
        virtual double intersect(Ray ray) = 0;
        virtual point getNormal(point p) = 0;
};

extern point checkerBoardPos;
extern point pos;
class Triangle {
    public:
    point a, b, c;
    Triangle(){
        a = point(0, 0, 0);
        b = point(0, 0, 0);
        c = point(0, 0, 0);
    }

    Triangle(point a, point b, point c){
        this->a = a;
        this->b = b;
        this->c = c;
    }

    void draw(){
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }

    point getNormal(){
        point ab = b - a;
        point ac = c - a;
        point normal = ab ^ ac;
        normal.normalize();
        return normal;
    }

    int pointInside(point p){
        point normal = getNormal();
        point ab = b - a;
        point bc = c - b;
        point ca = a - c;
        point ap = p - a;
        point bp = p - b;
        point cp = p - c;
        point cross1 = ab ^ ap;
        point cross2 = bc ^ bp;
        point cross3 = ca ^ cp;
        if ((cross1 * cross2) >= 0 && (cross2 * cross3) >= 0)
            return 1;
        return 0;
    }

    double intersect(Ray ray){
        point normal = getNormal();
        double t = -1 * ((ray.start - a) * normal) / (ray.dir * normal);
        if (t < 0)
            return -1;
        point p = ray.getPoint(t);
        point ab = b - a;
        point bc = c - b;
        point ca = a - c;
        point ap = p - a;
        point bp = p - b;
        point cp = p - c;
        point cross1 = ab ^ ap;
        point cross2 = bc ^ bp;
        point cross3 = ca ^ cp;
        if ((cross1 * cross2) >= 0 && (cross2 * cross3) >= 0)
            return t;
        return -1;
    }

    friend ostream& operator<<(ostream& os, const Triangle& triangle){
        os<<"Triangle: "<<endl;
        os<<"A: "<<triangle.a.x<<" "<<triangle.a.y<<" "<<triangle.a.z<<endl;
        os<<"B: "<<triangle.b.x<<" "<<triangle.b.y<<" "<<triangle.b.z<<endl;
        os<<"C: "<<triangle.c.x<<" "<<triangle.c.y<<" "<<triangle.c.z<<endl;
        return os;
    }
};

class CheckerBoard:public Object {
    public:
    double width;
    double leftX, leftY, rightX, rightY;

    CheckerBoard(double length, double coeff[3]){
        width = length;
        leftX = -length * 30;
        leftY = -length * 30;
        rightX = length * 30;
        rightY = length * 30;

        coEfficients[0] = coeff[0];
        coEfficients[1] = coeff[1];
        coEfficients[2] = 0;
        coEfficients[3] = coeff[2];
        shine = 0;
    }

    void draw() override{
        // cout<<width<<" "<<left<<endl;
        // if(checkerBoardPos.x - pos)
        double x_mov = (checkerBoardPos.x - pos.x) /(2*width);
        double y_mov = (checkerBoardPos.y - pos.y) /(2*width);
        leftX -= floor(x_mov) * width;
        rightX -= floor(x_mov) * width;
        leftY -= floor(y_mov) * width;
        rightY -= floor(y_mov) * width;
        checkerBoardPos.x -= floor(x_mov) * width;
        checkerBoardPos.y -= floor(y_mov) * width;
        // cout<<*this<<endl;

        for (double x = leftX; x < rightX; x+=width)
        {
            for (double y = leftY; y < rightY; y+=width)
            {
                // float x = i * width;
                // float y = j * width;
                int i = (int)((x) / width);
                int j = (int)((y) / width);

                glBegin(GL_QUADS);

                if ((i + j) % 2 == 0)
                {
                    glColor3f(0.5, 0.5, 0.5); // Light gray for even-indexed boxes
                }
                else
                {
                    glColor3f(0.25, 0.25, 0.25); // Dark gray for odd-indexed boxes
                }

                glVertex3f(x, y, 0);
                glVertex3f(x + width, y, 0);
                glVertex3f(x + width, y + width, 0);
                glVertex3f(x, y + width, 0);
                glEnd();
            }
        }
    }

    point getColor(point p) override{
        int x = (int)((p.x - leftX) / width);
        int y = (int)((p.y - leftY) / width);
        int z = x + y;
        //checkloadTexture
        if (loadTexture)
        {
            if (z % 2 == 0)
            {
                // return point(textureImageW.red_at(x, y) / 255.0, textureImageW.green_at(x, y) / 255.0, textureImageW.blue_at(x, y) / 255.0);
                double int_x = (p.x - leftX) - (x * width);
                double int_y = (p.y - leftY) - (y * width);
                int textureX = (int)(int_x * textureImageW.width() / width);
                int textureY = (int)(int_y * textureImageW.height() / width);
                //clamp textureX and textureY
                if (textureX < 0)
                    textureX = 0;
                if (textureX >= textureImageW.width())
                    textureX = textureImageW.width() - 1;
                if (textureY < 0)
                    textureY = 0;
                if (textureY >= textureImageW.height())
                    textureY = textureImageW.height() - 1;
                unsigned char r, g, b;
                textureImageW.get_pixel(textureX, textureY, r, g, b);
                return point(r / 255.0, g / 255.0, b / 255.0);    
            }
            else
            {
                double int_x = (p.x - leftX) - (x * width);
                double int_y = (p.y - leftY) - (y * width);
                int textureX = (int)(int_x * textureImageB.width() / width);
                int textureY = (int)(int_y * textureImageB.height() / width);
                //clamp textureX and textureY
                if (textureX < 0)
                    textureX = 0;
                if (textureX >= textureImageB.width())
                    textureX = textureImageB.width() - 1;
                if (textureY < 0)
                    textureY = 0;
                if (textureY >= textureImageB.height()) 
                    textureY = textureImageB.height() - 1;
                unsigned char r, g, b;
                textureImageB.get_pixel(textureX, textureY, r, g, b);
                return point(r / 255.0, g / 255.0, b / 255.0);
                // return point(textureImageB.red_at(x, y) / 255.0, textureImageB.green_at(x, y) / 255.0, textureImageB.blue_at(x, y) / 255.0);
            }
        }
        else
        {
            if (z % 2 == 0)
            {
                return point(1, 1, 1);
            }
            else
            {
                return point(0, 0, 0);
            }
        }
    }

    double intersect(Ray ray) override{
        double t = -1 * (ray.start.z / ray.dir.z);
        if (t < 0)
            return -1;
        point p = ray.getPoint(t);
        if (p.x < leftX || p.x > rightX || p.y < leftY || p.y > -leftY)
            return -1;
        return t;
    }

    point getNormal(point p) override{
        // return point(0, 0, 1);
        if (p.z > 0)
            return point(0, 0, 1);
        else
            return point(0, 0, -1);
    }

    friend ostream& operator<<(ostream& os, const CheckerBoard& checkerBoard){
        os<<"CheckerBoard: "<< endl;
        os<<"Width: "<<checkerBoard.width<<endl;
        os<<"Ambient  diffuse  specular  reflection: "<<checkerBoard.coEfficients[0]<<" "<<checkerBoard.coEfficients[1]<<" "<<checkerBoard.coEfficients[2]<<" "<<checkerBoard.coEfficients[3]<<endl;
    }
};

class Sphere : public Object
{
    public:
    point center;
    double radius;

    Sphere(point center, double radius,double color[3], double coeff[4], double shine)
    {
        this->center = center;
        this->radius = radius;
        coEfficients[0] = coeff[0];
        coEfficients[1] = coeff[1];
        coEfficients[2] = coeff[2];
        coEfficients[3] = coeff[3];
        this->color = point(color[0], color[1], color[2]);
        this->shine = shine;
    }

    void draw() override
    {
        glTranslatef(center.x, center.y, center.z);
        double r = color.x;
        double g = color.y;
        double b = color.z;
        struct point points[100][100];
        int i, j;
        double h, rd;
        int stacks = 20;
        int slices = 24;
        // generate points
        for (i = 0; i <= stacks; i++)
        {
            h = radius * sin(((double)i / (double)stacks) * (pi / 2));
            rd = radius * cos(((double)i / (double)stacks) * (pi / 2));
            for (j = 0; j <= slices; j++)
            {
                points[i][j].x = rd * cos(((double)j / (double)slices) * 2 * pi);
                points[i][j].y = rd * sin(((double)j / (double)slices) * 2 * pi);
                points[i][j].z = h;
            }
        }
        // cout << *this << endl;
        // draw quads using generated points
        for (i = 0; i < stacks; i++)
        {
            for (j = 0; j < slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    // upper hemisphere
                    glColor3f(r, g, b);
                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
                    // lower hemisphere
                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
                }
                glEnd();
            }
        }
        glTranslatef(-center.x, -center.y, -center.z);
    }

    point getColor(point p) override
    {
        return color;
    }

    double intersect(Ray ray) override
    {
        point start = ray.start;
        point dir = ray.dir;
        point oc = start - center;
        double a = dir * dir;
        double b = 2 * (oc * dir);
        double c = (oc * oc) - radius * radius;
        double d = b * b - 4 * a * c;
        if (d < 0)
            return -1;
        d = sqrt(d);
        double t1 = (-b + d) / (2 * a);
        double t2 = (-b - d) / (2 * a);
        if (t1 < 0 && t2 < 0)
            return -1;
        if (t1 < 0)
            return t2;
        if (t2 < 0)
            return t1;
        return min(t1, t2);
    }

    point getNormal(point p) override
    {
        // point normal = p - center;
        // normal.normalize();
        // return normal;

        point normal = p - center;
        normal.normalize();
        
        return normal;
    }

    friend ostream& operator<<(ostream& os, const Sphere& sphere){
        os<<"Sphere: "<<endl;
        os<<"Center: "<<sphere.center.x<<" "<<sphere.center.y<<" "<<sphere.center.z<<endl;
        os<<"Radius: "<<sphere.radius<<endl;
        os<<"Ambient  diffuse  specular  reflection: "<<sphere.coEfficients[0]<<" "<<sphere.coEfficients[1]<<" "<<sphere.coEfficients[2]<<" "<<sphere.coEfficients[3]<<endl;
        os<<"Color: "<<sphere.color.x<<" "<<sphere.color.y<<" "<<sphere.color.z<<endl;
        os<<"Shine: "<<sphere.shine<<endl;
        return os;
    }
};

class Cube : public Object {
    public:
    point lowerLeft;
    double side;
    point vertices[8];
    Triangle triangles[12];
    Cube(point lowerLeft, double side, double color[3], double coeff[4], double shine){
        this->lowerLeft = lowerLeft;
        this->side = side;
        coEfficients[0] = coeff[0];
        coEfficients[1] = coeff[1];
        coEfficients[2] = coeff[2];
        coEfficients[3] = coeff[3];
        this->color = point(color[0], color[1], color[2]);
        this->shine = shine;
        calculateVertices();
        findTriangles();
    }

    void calculateVertices(){
        vertices[0] = lowerLeft;
        vertices[1] = lowerLeft + point(side, 0, 0);
        vertices[2] = lowerLeft + point(side, side, 0);
        vertices[3] = lowerLeft + point(0, side, 0);
        vertices[4] = lowerLeft + point(0, 0, side);
        vertices[5] = lowerLeft + point(side, 0, side);
        vertices[6] = lowerLeft + point(side, side, side);
        vertices[7] = lowerLeft + point(0, side, side);
    }

    void findTriangles(){
        triangles[0] = Triangle(vertices[0], vertices[1], vertices[2]);
        triangles[1] = Triangle(vertices[0], vertices[2], vertices[3]);
        triangles[2] = Triangle(vertices[0], vertices[4], vertices[5]);
        triangles[3] = Triangle(vertices[0], vertices[5], vertices[1]);
        triangles[4] = Triangle(vertices[1], vertices[5], vertices[6]);
        triangles[5] = Triangle(vertices[1], vertices[6], vertices[2]);
        triangles[6] = Triangle(vertices[2], vertices[6], vertices[7]);
        triangles[7] = Triangle(vertices[2], vertices[7], vertices[3]);
        triangles[8] = Triangle(vertices[3], vertices[7], vertices[4]);
        triangles[9] = Triangle(vertices[3], vertices[4], vertices[0]);
        triangles[10] = Triangle(vertices[4], vertices[7], vertices[6]);
        triangles[11] = Triangle(vertices[4], vertices[6], vertices[5]);
    }

    void draw() override{
        // cout << *this << endl;
        glColor3f(color.x, color.y, color.z);
        for (int i = 0; i < 12; i++)
            triangles[i].draw();
        // glPushMatrix();
        // {
        //     glTranslatef(lowerLeft.x, lowerLeft.y, lowerLeft.z);
        //     glScalef(side, side, side);

        //     double r = color.x;
        //     double g = color.y;
        //     double b = color.z;
        //     glBegin(GL_QUADS); // Begin drawing the color cube with 6 quads
        //     // Top face (y = 1.0f)
        //     // Define vertices in counter-clockwise (CCW) order with normal pointing out
        //     glColor3f(r, g, b);

        //     // First Quad
        //     glBegin(GL_TRIANGLES);
        //     glVertex3f(1.0f, 1.0f, -1.0f);
        //     glVertex3f(-1.0f, 1.0f, -1.0f);
        //     glVertex3f(-1.0f, 1.0f, 1.0f);

        //     glVertex3f(-1.0f, 1.0f, 1.0f);
        //     glVertex3f(1.0f, 1.0f, 1.0f);
        //     glVertex3f(1.0f, 1.0f, -1.0f);
        //     glEnd();

        //     // Second Quad
        //     glBegin(GL_TRIANGLES);
        //     glVertex3f(1.0f, -1.0f, 1.0f);
        //     glVertex3f(-1.0f, -1.0f, 1.0f);
        //     glVertex3f(-1.0f, -1.0f, -1.0f);

        //     glVertex3f(-1.0f, -1.0f, -1.0f);
        //     glVertex3f(1.0f, -1.0f, -1.0f);
        //     glVertex3f(1.0f, -1.0f, 1.0f);
        //     glEnd();

        //     // Third Quad
        //     glBegin(GL_TRIANGLES);
        //     glVertex3f(1.0f, 1.0f, 1.0f);
        //     glVertex3f(-1.0f, 1.0f, 1.0f);
        //     glVertex3f(-1.0f, -1.0f, 1.0f);

        //     glVertex3f(-1.0f, -1.0f, 1.0f);
        //     glVertex3f(1.0f, -1.0f, 1.0f);
        //     glVertex3f(1.0f, 1.0f, 1.0f);
        //     glEnd();

        //     // Fourth Quad
        //     glBegin(GL_TRIANGLES);
        //     glVertex3f(1.0f, -1.0f, -1.0f);
        //     glVertex3f(-1.0f, -1.0f, -1.0f);
        //     glVertex3f(-1.0f, 1.0f, -1.0f);

        //     glVertex3f(-1.0f, 1.0f, -1.0f);
        //     glVertex3f(1.0f, 1.0f, -1.0f);
        //     glVertex3f(1.0f, -1.0f, -1.0f);
        //     glEnd();

        //     // Fifth Quad
        //     glBegin(GL_TRIANGLES);
        //     glVertex3f(-1.0f, 1.0f, 1.0f);
        //     glVertex3f(-1.0f, 1.0f, -1.0f);
        //     glVertex3f(-1.0f, -1.0f, -1.0f);

        //     glVertex3f(-1.0f, -1.0f, -1.0f);
        //     glVertex3f(-1.0f, -1.0f, 1.0f);
        //     glVertex3f(-1.0f, 1.0f, 1.0f);
        //     glEnd();

        //     // Sixth Quad
        //     glBegin(GL_TRIANGLES);
        //     glVertex3f(1.0f, 1.0f, -1.0f);
        //     glVertex3f(1.0f, 1.0f, 1.0f);
        //     glVertex3f(1.0f, -1.0f, 1.0f);

        //     glVertex3f(1.0f, -1.0f, 1.0f);
        //     glVertex3f(1.0f, -1.0f, -1.0f);
        //     glVertex3f(1.0f, 1.0f, -1.0f);
        //     glEnd();

        //     glEnd();
        // }
        // glPopMatrix();
    }
    
    point getColor(point p) override{
        return color;
    }

    double intersect(Ray ray) override{
        vector<double> vc_t;
        for (int i = 0; i < 12; i++)
        {
            double t = triangles[i].intersect(ray);
            if (t > 0)
                vc_t.push_back(t);
        }
        if (vc_t.size() == 0)
            return -1;
        sort(vc_t.begin(), vc_t.end());
        return vc_t[0];
    }

    point getNormal(point p) override{
        for (int i = 0; i < 12; i++)
        {
            if (triangles[i].pointInside(p))
                return triangles[i].getNormal();
        }
        return point(0, 0, 0);
    }

    friend ostream& operator<<(ostream& os, const Cube& cube){
        os<<"Cube: "<<endl;
        os<<"LowerLeft: "<<cube.lowerLeft.x<<" "<<cube.lowerLeft.y<<" "<<cube.lowerLeft.z<<endl;
        os<<"Side: "<<cube.side<<endl;
        os<<"Ambient  diffuse  specular  reflection: "<<cube.coEfficients[0]<<" "<<cube.coEfficients[1]<<" "<<cube.coEfficients[2]<<" "<<cube.coEfficients[3]<<endl;
        os<<"Color: "<<cube.color.x<<" "<<cube.color.y<<" "<<cube.color.z<<endl;
        os<<"Shine: "<<cube.shine<<endl;
        return os;
    }
};

class Pyramid : public Object {
    public:
    point lowerLeft;
    double width, height;
    point vertices[5];
    Triangle triangles[6];

    Pyramid(point lowerLeft, double width, double height, double color[3], double coeff[4], double shine){
        this->lowerLeft = lowerLeft;
        this->width = width;
        this->height = height;
        coEfficients[0] = coeff[0];
        coEfficients[1] = coeff[1];
        coEfficients[2] = coeff[2];
        coEfficients[3] = coeff[3];
        this->color = point(color[0], color[1], color[2]);
        this->shine = shine;
        calculateVertices();
        findTriangles();
        // for (int i = 0; i < 6; i++)
        //     cout << triangles[i] << endl;
    }

    void draw() override{
        // cout << *this << endl;
        glColor3f(color.x, color.y, color.z);
        for (int i = 0; i < 6; i++)
            triangles[i].draw();
        // glPushMatrix();
        // {
        //     glTranslatef(lowerLeft.x, lowerLeft.y, lowerLeft.z);
        //     glScalef(width, width, height);

        //     double r = color.x;
        //     double g = color.y;
        //     double b = color.z;
        //     glBegin(GL_TRIANGLES);
        //     // Triangle 1
        //     glColor3f(r, g, b);
        //     glVertex3f(0.0f, 0.0f, 0.0f);
        //     glVertex3f(1.0f, 0.0f, 0.0f);
        //     glVertex3f(1.0f, 1.0f, 0.0f);

        //     // Triangle 2
        //     glVertex3f(1.0f, 1.0f, 0.0f);
        //     glVertex3f(0.0f, 1.0f, 0.0f);
        //     glVertex3f(0.0f, 0.0f, 0.0f);

        //     glVertex3f(0.0f, 0.0f, 0.0f);
        //     glVertex3f(1.0f, 0.0f, 0.0f);
        //     glVertex3f(0.5f, 0.5f, 1.0f);

        //     glVertex3f(0.0f, 0.0f, 0.0f);
        //     glVertex3f(0.5f, 0.5f, 1.0f);
        //     glVertex3f(0.0f, 1.0f, 0.0f);

        //     glVertex3f(0.0f, 1.0f, 0.0f);
        //     glVertex3f(0.5f, 0.5f, 1.0f);
        //     glVertex3f(1.0f, 1.0f, 0.0f);

        //     glVertex3f(1.0f, 1.0f, 0.0f);
        //     glVertex3f(0.5f, 0.5f, 1.0f);
        //     glVertex3f(1.0f, 0.0f, 0.0f);
        //     glEnd(); // Done drawing the pyramid
        // }
        // glPopMatrix();
    }

    point getColor(point p) override
    {
        return color;
    }

    double intersect(Ray ray) override
    {
        vector<double> vc_t;
        for (int i = 0; i < 6; i++)
        {
            double t = triangles[i].intersect(ray);
            if (t > 0)
                vc_t.push_back(t);
        }
        if (vc_t.size() == 0)
            return -1;
        double t = vc_t[0];
        for (int i = 1; i < vc_t.size(); i++)
        {
            if (vc_t[i] < t)
                t = vc_t[i];
        }
        return t;
    }

    point getNormal(point p) override
    {
        for (int i = 0; i < 6; i++)
        {
            if (triangles[i].pointInside(p))
                return triangles[i].getNormal();
        }
    }

    void calculateVertices(){
        vertices[0] = lowerLeft;
        vertices[1] = lowerLeft + point(width, 0, 0);
        vertices[2] = lowerLeft + point(width, width, 0);
        vertices[3] = lowerLeft + point(0, width, 0);
        vertices[4] = lowerLeft + point(width/2, width/2, height);
    }
    
    void findTriangles(){
        triangles[0] = Triangle(vertices[0], vertices[1], vertices[4]);
        triangles[1] = Triangle(vertices[1], vertices[2], vertices[4]);
        triangles[2] = Triangle(vertices[2], vertices[3], vertices[4]);
        triangles[3] = Triangle(vertices[3], vertices[0], vertices[4]);
        triangles[4] = Triangle(vertices[0], vertices[1], vertices[2]);
        triangles[5] = Triangle(vertices[0], vertices[2], vertices[3]);
    }

    friend ostream& operator<<(ostream& os, const Pyramid& pyramid){
        os<<"Pyramid: "<<endl;
        os<<"LowerLeft: "<<pyramid.lowerLeft.x<<" "<<pyramid.lowerLeft.y<<" "<<pyramid.lowerLeft.z<<endl;
        os<<"Width: "<<pyramid.width<<endl;
        os<<"Height: "<<pyramid.height<<endl;
        os<<"Ambient  diffuse  specular  reflection: "<<pyramid.coEfficients[0]<<" "<<pyramid.coEfficients[1]<<" "<<pyramid.coEfficients[2]<<" "<<pyramid.coEfficients[3]<<endl;
        os<<"Color: "<<pyramid.color.x<<" "<<pyramid.color.y<<" "<<pyramid.color.z<<endl;
        os<<"Shine: "<<pyramid.shine<<endl;
        return os;
    }
};

class Light{
    public:
    point pos;
    double falloff;

    Light(point pos, double falloff){
        this->pos = pos;
        this->falloff = falloff;
    }

    friend ostream& operator<<(ostream& os, const Light& light){
        os<<"Light: "<<light.pos.x<<" "<<light.pos.y<<" "<<light.pos.z<<" "<<light.falloff<<endl;
        return os;
    }
};

class SpotLight{
    public:
    point pos;
    double falloff;
    point dir;
    double cutoff;

    SpotLight(point pos, double falloff, point dir, double cutoff){
        this->pos = pos;
        this->falloff = falloff;
        this->dir = dir;
        this->cutoff = cutoff;
    }

    friend ostream& operator<<(ostream& os, const SpotLight& light){
        os<<"SpotLight: "<<light.pos.x<<" "<<light.pos.y<<" "<<light.pos.z<<" "<<light.falloff<<" "<<light.dir.x<<" "<<light.dir.y<<" "<<light.dir.z<<" "<<light.cutoff<<endl;
        return os;
    }
};