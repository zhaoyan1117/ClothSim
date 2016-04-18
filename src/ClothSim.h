#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#ifdef OSX
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <omp.h>

#define PI 3.14159265
inline float sqr(float x) { return x*x; }

#define DELTA_T 0.005

using namespace std;

class Particle {
  public:
    bool movable;
    float mass;
    glm::vec2 planeCoord;
    glm::vec3 pos;
    glm::vec3 prev_pos;
    glm::vec3 vel;
    glm::vec3 acceleration;
    glm::vec3 normal;
    glm::vec3 color;

    //For .obj file output
    int index;

    Particle(glm::vec3 p, glm::vec2 pc, int i) {
        pos = p;
        prev_pos = p;
        planeCoord = pc;

        movable = true;
        mass = 1.0f;
        vel = glm::vec3(0.0f, 0.0f, 0.0f);
        acceleration = glm::vec3(0.0f, 0.0f, 0.0f);
        normal = glm::vec3(0.0f, 0.0f, 0.0f);
        color = glm::vec3();

        index = i;
    }

    Particle() {}

    void addForce(glm::vec3 force) {
        // Sicen mass is 1.0f;
        acceleration += force;

        // Used when mass is not 1.0f;
        //acceleration += (force/mass);
    }

    void addNormal(glm::vec3 n) {
        normal += glm::normalize(n);
    }

    void resetNormal() {
        normal = glm::vec3(0.0f, 0.0f, 0.0f);
    }

    void move(glm::vec3 dir) {
        pos = pos + dir;
    }

    void timeStep() {
        if (movable) {

            //Verlet Intergration.
            glm::vec3 temp = pos;
            pos += (pos-prev_pos) + acceleration*(float)DELTA_T;
            prev_pos = temp;
            vel = pos-prev_pos;
            
            //Euler Integration.
            //vel += acceleration * (float) DELTA_T;
            //pos += vel * (float) DELTA_T;

            //cout << "acceleration: " << acceleration.x << " " << acceleration.y << " " << acceleration.z << endl;
            //cout << "pos: " << pos.x << " " << pos.y << " " << pos.z << endl;

            color = glm::vec3((acceleration.x)*30.0f, (acceleration.y)*30.0f, (acceleration.z)*30.0f);
            acceleration = glm::vec3(0.0f, 0.0f, 0.0f);
        }
    }
};

class TriangleMesh {
  public:
    Particle *p1, *p2, *p3;
    float kst, ksh, dst;
    float initW, initH;
    glm::mat2 theInverse;

    float A;

    glm::vec3 uix, uiy, uiz;
    glm::vec3 ujx, ujy, ujz;
    glm::vec3 ukx, uky, ukz;

    glm::vec3 vix, viy, viz;
    glm::vec3 vjx, vjy, vjz;
    glm::vec3 vkx, vky, vkz;

    TriangleMesh() {}

    TriangleMesh(Particle* pp1, Particle* pp2, Particle* pp3, float st, float sh, float ds, float w, float h) {
        p1 = pp1; p2 = pp2; p3 = pp3;
        kst = st; ksh = sh; dst = ds;
        initW = w; initH = h;

        //Setup.
        float deltaU1 = p2->planeCoord.x - p1->planeCoord.x;
        float deltaV1 = p2->planeCoord.y - p1->planeCoord.y;
        float deltaU2 = p3->planeCoord.x - p1->planeCoord.x;
        float deltaV2 = p3->planeCoord.y - p1->planeCoord.y;

        glm::mat2 temp(deltaU1, deltaV1, deltaU2, deltaV2);
        theInverse = glm::inverse(temp);

        float u1v2u2v1 = deltaU1*deltaV2 - deltaU2*deltaV1;
        A = 0.5f * (u1v2u2v1);

        float ui, uj, uk, vi, vj, vk;
        if (u1v2u2v1 != 0.0f) {
            ui = (deltaV1 - deltaV2) / (u1v2u2v1);
            uj = (deltaV2) / (u1v2u2v1);
            uk = (-deltaV1) / (u1v2u2v1);
            vi = (deltaU2 - deltaU1) / (u1v2u2v1);
            vj = (-deltaU2) / (u1v2u2v1);
            vk = (deltaU1) / (u1v2u2v1);
        } else {
            ui = 0.0f;
            uj = 0.0f;
            uk = 0.0f;
            vi = 0.0f;
            vj = 0.0f;
            vk = 0.0f;
        }

        uix = glm::vec3(ui, 0.0f, 0.0f); uiy = glm::vec3(0.0f, ui, 0.0f); uiz = glm::vec3(0.0f, 0.0f, ui);
        ujx = glm::vec3(uj, 0.0f, 0.0f); ujy = glm::vec3(0.0f, uj, 0.0f); ujz = glm::vec3(0.0f, 0.0f, uj);
        ukx = glm::vec3(uk, 0.0f, 0.0f); uky = glm::vec3(0.0f, uk, 0.0f); ukz = glm::vec3(0.0f, 0.0f, uk);

        vix = glm::vec3(vi, 0.0f, 0.0f); viy = glm::vec3(0.0f, vi, 0.0f); viz = glm::vec3(0.0f, 0.0f, vi);
        vjx = glm::vec3(vj, 0.0f, 0.0f); vjy = glm::vec3(0.0f, vj, 0.0f); vjz = glm::vec3(0.0f, 0.0f, vj);
        vkx = glm::vec3(vk, 0.0f, 0.0f); vky = glm::vec3(0.0f, vk, 0.0f); vkz = glm::vec3(0.0f, 0.0f, vk);
    }

    void timeStep() {
        glm::vec3 deltaX1 = p2->pos - p1->pos;
        glm::vec3 deltaX2 = p3->pos - p1->pos;

        glm::mat2x3 deltaX12(1);
        deltaX12 = glm::column(deltaX12, 0, deltaX1);
        deltaX12 = glm::column(deltaX12, 1, deltaX2);

        glm::mat2x3 wUV = deltaX12 * theInverse;
        glm::vec3 wU = glm::column(wUV, 0);
        glm::vec3 wV = glm::column(wUV, 1);
        float wUl = glm::length(wU); float wVl = glm::length(wV);

        glm::vec2 Cst = A * glm::vec2(wUl-initW, wVl-initH);

        float factorU = (-A/wUl);
        float factorV = (-A/wVl);

        // Stretch(Damping) Force for particle p1(i).
        glm::vec3 tempU = glm::vec3(glm::dot(uix, wU),
                                    glm::dot(uiy, wU),
                                    glm::dot(uiz, wU));
        tempU *= factorU;

        glm::vec3 tempV = glm::vec3(glm::dot(vix, wV),
                                    glm::dot(viy, wV),
                                    glm::dot(viz, wV));
        tempV *= factorV;

        glm::mat2x3 forceST(1);
        forceST = glm::column(forceST, 0, tempU);
        forceST = glm::column(forceST, 1, tempV);

        glm::vec3 fSTi = forceST * Cst;
        p1->addForce(fSTi*kst);

        glm::mat3x2 dampST(1);
        dampST = glm::row(dampST, 0, tempU);
        dampST = glm::row(dampST, 1, tempV);
        glm::vec2 dampST2 = dampST * p1->vel;
        glm::vec3 dSTi = forceST*dampST2;
        p1->addForce(dSTi*(-dst));


        // Stretch(Damping) Force for particle p2(j).
        tempU = glm::vec3(glm::dot(ujx, wU),
                          glm::dot(ujy, wU),
                          glm::dot(ujz, wU));
        tempU *= factorU;

        tempV = glm::vec3(glm::dot(vjx, wV),
                          glm::dot(vjy, wV),
                          glm::dot(vjz, wV));
        tempV *= factorV;

        forceST = glm::column(forceST, 0, tempU);
        forceST = glm::column(forceST, 1, tempV);

        glm::vec3 fSTj = forceST * Cst;
        p2->addForce(fSTj*kst);

        dampST = glm::row(dampST, 0, tempU);
        dampST = glm::row(dampST, 1, tempV);
        dampST2 = dampST * p2->vel;
        glm::vec3 dSTj = forceST*dampST2;
        p2->addForce(dSTj*(-dst));

        // Stretch(Damping) Force for particle p3(k).
        tempU = glm::vec3(glm::dot(ukx, wU),
                          glm::dot(uky, wU),
                          glm::dot(ukz, wU));
        tempU *= factorU;

        tempV = glm::vec3(glm::dot(vkx, wV),
                          glm::dot(vky, wV),
                          glm::dot(vkz, wV));
        tempV *= factorV;

        forceST = glm::column(forceST, 0, tempU);
        forceST = glm::column(forceST, 1, tempV);

        glm::vec3 fSTk = forceST * Cst;
        p3->addForce(fSTk*kst);

        dampST = glm::row(dampST, 0, tempU);
        dampST = glm::row(dampST, 1, tempV);
        dampST2 = dampST * p3->vel;
        glm::vec3 dSTk = forceST*dampST2;
        p3->addForce(dSTk*(-dst));

        // For shear force calculation.
        float Csh = A * glm::dot(wV, wU);
        float shearFactor = -ksh*A*Csh;
        // For p1(i).
        glm::vec3 shearForce = glm::vec3(glm::dot(uix,wV)+glm::dot(vix,wU),
                                         glm::dot(uiy,wV)+glm::dot(viy,wU),
                                         glm::dot(uiz,wV)+glm::dot(viz,wU));
        shearForce *= shearFactor;
        p1->addForce(shearForce);

        // For p2(j).
        shearForce = glm::vec3(glm::dot(ujx,wV)+glm::dot(vjx,wU),
                               glm::dot(ujy,wV)+glm::dot(vjy,wU),
                               glm::dot(ujz,wV)+glm::dot(vjz,wU));
        shearForce *= shearFactor;
        p2->addForce(shearForce);

        // For p3(k).
        shearForce = glm::vec3(glm::dot(ukx,wV)+glm::dot(vkx,wU),
                               glm::dot(uky,wV)+glm::dot(vky,wU),
                               glm::dot(ukz,wV)+glm::dot(vkz,wU));
        shearForce *= shearFactor;
        p3->addForce(shearForce);
    }
};

class EdgeSpring {
  public:
    Particle *p1, *p2, *p3, *p4;
    float ke, kd;

    EdgeSpring(Particle *pp1, Particle *pp2, Particle *pp3, Particle *pp4, float e, float d) {
        p1 = pp1, p2 = pp2, p3 = pp3, p4 = pp4;
        ke = e; kd = d;
    }

    EdgeSpring() {}

    void timeStep() {
        glm::vec3 e13 = p1->pos - p3->pos;
        glm::vec3 e14 = p1->pos - p4->pos;
        glm::vec3 e23 = p2->pos - p3->pos;
        glm::vec3 e24 = p2->pos - p4->pos;
        glm::vec3 E = p4->pos - p3->pos;

        glm::vec3 N1 = glm::cross(e13, e14);
        glm::vec3 N2 = glm::cross(e24, e23);

        float El = glm::length(E);
        float N1l2 = sqr(glm::length(N1));
        float N2l2 = sqr(glm::length(N2));
        
        glm::vec3 u1 = N1 * (El/N1l2);
        glm::vec3 u2 = N2 * (El/N2l2);
        glm::vec3 u3 = N1 * glm::dot(e14, E)/(El*N1l2) + N2 * glm::dot(e24, E)/(El*N2l2);
        glm::vec3 u4 = - N1 * glm::dot(e13, E)/(El*N1l2) - N2 * glm::dot(e23, E)/(El*N2l2);

        glm::vec3 n1 = glm::normalize(N1); glm::vec3 n2 = glm::normalize(N2);
        float sinHalf = sqrt((1.0f-glm::dot(n1, n2))/2.0f);
        if (sinHalf != 0.0f) {
            if (glm::dot(glm::cross(n1, n2), glm::normalize(E)) < 0) { sinHalf = -sinHalf; }
            float bendFactor = ke*sqr(El)*sinHalf/(N1l2+N2l2);

            if (bendFactor != bendFactor) {bendFactor = 0.0f;}

            p1->addForce(u1*bendFactor);
            p2->addForce(u2*bendFactor);
            p3->addForce(u3*bendFactor);
            p4->addForce(u4*bendFactor);
        }

        // Damping forces.
        float dThetaDT = glm::dot(u1, p1->vel) + glm::dot(u2, p2->vel) + glm::dot(u3, p3->vel) + glm::dot(u4, p4->vel);
        float bdFactor = -kd*El*dThetaDT;
        if (bdFactor != bdFactor) {bdFactor = 0.0f;}
        p1->addForce(u1*bdFactor);
        p2->addForce(u2*bdFactor);
        p3->addForce(u3*bdFactor);
        p4->addForce(u4*bdFactor);

    }
};

class Cloth {
  private:
    std::vector<Particle> particlesBACKUP;

  public:
    int widthNum, heightNum;
    float initW, initH;
    std::vector<Particle> particles;

    std::vector<EdgeSpring> ess;
    std::vector<TriangleMesh> tms;

    float ksh, kst, dst;
    float ke, kd;

    // For aerodynamic.
    glm::vec3 windVel;

    // For ball collision.
    bool ballCol;
    glm::vec3 center;
    float radius;
    glm::vec3 speed;

    // For cube collision.
    bool cubeCol;
    glm::vec3 cubeMax, cubeMin;

    // For floor collision.
    bool hasFloor;
    float floor;

    // For .obj file output.
    int numOfParticles;

    Cloth() {}

    void makeEdgeSpring (Particle *p1, Particle *p2, Particle *p3, Particle *p4) {
        ess.push_back(EdgeSpring(p1, p2, p3, p4, ke, kd));
    }

    void makeTriangleMesh (Particle *p1, Particle *p2, Particle *p3) {
        tms.push_back(TriangleMesh(p1, p2, p3, kst, ksh, dst, initW, initH));
    }

    Particle* getParticle(int x, int y) {
        return &(particles[y*widthNum + x]);
    }

    void setMovable(int x, int y, bool m) {
        getParticle(x, y)->movable = m; 
    }

    void drop() {
        for(int x=0; x<widthNum; x++) {
            for(int y=0; y<heightNum; y++) {
                setMovable(x, y, true);
            }
        }
    }
    
    void pull() {
        
        float midW = widthNum/2.0f;
        float midH = heightNum/2.0f;
        float wExt = 0.2f*initW/midW;
        float hExt = 0.2f*initH/midH;
        float mExt = (wExt + hExt)/2.0f;

        for(int x=0; x<widthNum; x++) {
            for(int y=0; y<heightNum; y++) {
                getParticle(x, y)->move(glm::vec3( abs((float)x-midW)*wExt, ((abs((float)x-midW)+abs((float)y-midH))/2.0f)*mExt, abs((float)y-midH)*hExt));
            }
        }

    }

    Cloth(float width, float height, int wn, int hn, float st, float sh, float ds, float e, float d, bool hang, float locHeight) {
        kst = st; ksh = sh; dst = ds;
        ke = e; kd = d;
        widthNum = wn; heightNum = hn;
        particles.resize(widthNum*heightNum);
        windVel = glm::vec3();
        hasFloor = false;

        // Initialize ball collision information.
        ballCol = false;
        center = glm::vec3();
        radius = 0.0f;

        // Initialize cube collision information.
        cubeCol = false;
        cubeMin = glm::vec3(); cubeMax = glm::vec3();

        numOfParticles = widthNum*heightNum;
        int curIndex = 1;

        initW = width / wn;
        initH = height / hn;
        // Creating Particles.
        if (hang) {
            for(int x=0; x<widthNum; x++) {
                for(int y=0; y<heightNum; y++) {
                    glm::vec3 pos = glm::vec3((float)x*initW,
                                    locHeight-(float)y*initH,
                                    0.0f);
                    glm::vec2 pcPos = glm::vec2(x,y);
                    particles[y*widthNum+x]= Particle(pos, pcPos, curIndex);
                    curIndex+=1;
                }
            }
        } else {
            for(int x=0; x<widthNum; x++) {
                for(int y=0; y<heightNum; y++) {
                    glm::vec3 pos = glm::vec3((float)x*initW,
                                    locHeight,
                                    (float)y*initH);
                    glm::vec2 pcPos = glm::vec2(x, y);
                    particles[y*widthNum+x]= Particle(pos, pcPos, curIndex);
                    curIndex+=1;
                }
            }
        }

        for(int x=0; x<widthNum-1; x++) {
            for(int y=0; y<heightNum-1; y++) {
                if ((x+y)%2 == 0){
                    makeTriangleMesh(getParticle(x,y), getParticle(x,y+1), getParticle(x+1,y));
                    makeTriangleMesh(getParticle(x+1,y+1), getParticle(x+1,y), getParticle(x,y+1));
                    makeEdgeSpring(getParticle(x,y), getParticle(x+1,y+1), getParticle(x,y+1), getParticle(x+1,y));
                } else {
                    makeTriangleMesh(getParticle(x,y), getParticle(x,y+1), getParticle(x+1,y+1));
                    makeTriangleMesh(getParticle(x,y), getParticle(x+1,y+1), getParticle(x+1,y));
                    makeEdgeSpring(getParticle(x,y+1), getParticle(x+1,y), getParticle(x,y), getParticle(x+1,y+1));
                }
            }
        }

        for(int x=0; x<widthNum-1; x++) {
            for(int y=0; y<heightNum-1; y++) {
                if (x<widthNum-2) {
                    if ((x+y)%2 == 0) {
                        makeEdgeSpring(getParticle(x,y+1), getParticle(x+2,y+1), getParticle(x+1,y), getParticle(x+1,y+1));
                    } else {
                        makeEdgeSpring(getParticle(x,y), getParticle(x+2,y), getParticle(x+1,y), getParticle(x+1,y+1));
                    }
                }
                if (y<heightNum-2) {
                    if ((x+y)%2 == 0) {
                        makeEdgeSpring(getParticle(x+1,y), getParticle(x+1,y+2), getParticle(x,y+1), getParticle(x+1,y+1));
                    } else {
                        makeEdgeSpring(getParticle(x,y), getParticle(x,y+2), getParticle(x,y+1), getParticle(x+1,y+1));
                    }
                }
            }
        }

        // Making unmovable points.
        if (hang) {
            for(int i=0;i<3; i++)
            {
                getParticle(0+i ,0)->movable = false; 
                getParticle(widthNum-1-i ,0)->movable = false;
                getParticle(widthNum/2-3+i, 0)->movable = false;
            }
        }
        
        // Save backups.
        particlesBACKUP = particles;

    }

    // For intersection tests.
    bool triangleIntersect(glm::vec3 old, glm::vec3 cur, float& intersect, glm::vec3 A, glm::vec3 B, glm::vec3 C) {
        glm::vec3 D = cur - old;
        float a, b, c, d, e, f, g, h, i, j, k, l, M, t, gamma, beta;
        float tMin = 0.0f, tMax = 1.0f;

        a = A.x - B.x;
        b = A.y - B.y;
        c = A.z - B.z;
        d = A.x - C.x;
        e = A.y - C.y;
        f = A.z - C.z;
        g = D.x;
        h = D.y;
        i = D.z;
        j = A.x - old.x;
        k = A.y - old.y;
        l = A.z - old.z;
        M = a*(e*i - h*f) + b*(g*f - d*i) + c*(d*h - e*g);

        t = (f*(a*k - j*b) + e*(j*c - a*l) + d*(b*l - k*c))/(-M);
        if (t < tMin or t > tMax) {
            return false;
        }

        gamma = (i*(a*k - j*b) + h*(j*c - a*l) + g*(b*l - k*c))/M;
        if (gamma < 0 or gamma > 1) {
            return false;
        }
        beta = (j*(e*i - h*f) + k*(g*f - d*i) + l*(d*h - e*g))/M;
        if (beta < 0 or beta > (1-gamma)) {
            return false;
        }

        // Check NaN.
        if (t != t) {
            t = 0.0f;
        }

        intersect = t;
        return true;
    }

    bool ballIntersect(glm::vec3 old, glm::vec3 cur, float& intersect) {
        glm::vec3 dir = cur - old;

        glm::vec3 eMinusC = old - center;
        float B_2 = sqr(glm::dot(dir, eMinusC));
        float AC_4 = glm::dot(dir, dir) * (glm::dot(eMinusC, eMinusC) - sqr(radius));
        float dis = B_2-AC_4;

        // Tangent and cross are both intersect.
        if (dis < 0) { return false; }

        float minusB = -glm::dot(dir, eMinusC);
        float D_2 = glm::dot(dir, dir);

        float t;
        t = (minusB-sqrt(dis)) / D_2;
        
        // Check for Nan.
        if (t != t) {t = 0.0f;}

        if (!(0.0f <= t and t <= 1.0f)) { return false; }

        intersect = t;

        return true;
    }

    bool hitSurface(glm::vec3 old, glm::vec3 cur, float& intersect, glm::vec3 v1, glm::vec3 v2, glm::vec3 v3, glm::vec3 v4) {
        float v123Intersect, v134Intersect;

        bool hit123 = triangleIntersect(old, cur, v123Intersect, v1, v2, v3);
        bool hit134 = triangleIntersect(old, cur, v134Intersect, v1, v3, v4);

        if (hit123) {
            intersect = v123Intersect;
            return hit123;
        } else if (hit134) {
            intersect = v134Intersect;
            return hit134;
        } else {
            return false;
        }
    }

    // For floor.
    bool hitFloor(glm::vec3 old, glm::vec3 cur, float& intersect) {
        glm::vec3 A = glm::vec3(-1000.0f, floor, -1000.0f);
        glm::vec3 B = glm::vec3(-1000.0f, floor, 1000.0f);
        glm::vec3 C = glm::vec3(1000.0f, floor, 1000.0f);
        glm::vec3 D = glm::vec3(1000.0f, floor, -1000.0f);
        return hitSurface(old, cur, intersect, A, B, C, D);
    }

    void setFloor(float f) {
        hasFloor = true;
        floor = f;
    }

    // For gravity.
    void addForce(glm::vec3 force) {
        for(std::vector<Particle>::size_type i = 0; i < particles.size(); i++) {
            particles[i].addForce(force);
        }
    }

    // For aerodynamic.
    glm::vec3 getNormal(Particle *p1,Particle *p2,Particle *p3) {
        glm::vec3 v1 = p2->pos - p1->pos;
        glm::vec3 v2 = p3->pos - p1->pos;

        return glm::normalize(glm::cross(v1, v2));
    }

    float getArea(Particle *p1,Particle *p2,Particle *p3) {
        glm::vec3 v1 = p2->pos - p1->pos;
        glm::vec3 v2 = p3->pos - p1->pos;

        return 0.5f * glm::cross(v1, v2).length();
    }

    void addWindForce(Particle *p1,Particle *p2,Particle *p3) {
        float density = 0.5f;
        float cod = 1.0f;

        glm::vec3 normal = getNormal(p1,p2,p3);

        glm::vec3 vSurface = (p1->vel + p2->vel + p3->vel) / 3.0f;
        glm::vec3 vRelative = vSurface - windVel;

        float a0 = getArea(p1, p2, p3);
        float a = (a0/vRelative.length())* glm::dot(vRelative, normal);

        glm::vec3 force = normal * (-0.5f * density * cod * a * sqr(vRelative.length()));
        p1->addForce(force);
        p2->addForce(force);
        p3->addForce(force);
    }

    void windForce() {
        for(int x = 0; x<widthNum - 1; x++) {
            for(int y=0; y<heightNum - 1; y++) {
                addWindForce(getParticle(x+1,y), getParticle(x,y), getParticle(x,y+1));
                addWindForce(getParticle(x+1,y+1), getParticle(x+1,y), getParticle(x,y+1));
            }
        }
        windVel = glm::vec3();
    }

    void updateWindForce(glm::vec3 wV) {
        windVel = wV;
    }

    // For drawing.
    void drawTriangle(Particle *p1, Particle *p2, Particle *p3, glm::vec3 color) {
        glColor3f(color.x, color.y, color.z);
        
        glm::vec3 p1n = glm::normalize(p1->normal);
        glNormal3f(p1n.x, p1n.y, p1n.z);
        glVertex3f(p1->pos.x, p1->pos.y, p1->pos.z);

        glm::vec3 p2n = glm::normalize(p2->normal);
        glNormal3f(p2n.x, p2n.y, p2n.z);
        glVertex3f(p2->pos.x, p2->pos.y, p2->pos.z);

        glm::vec3 p3n = glm::normalize(p3->normal);
        glNormal3f(p3n.x, p3n.y, p3n.z);
        glVertex3f(p3->pos.x, p3->pos.y, p3->pos.z);
    }

    void draw() {
        #pragma omp parallel for
        for(std::vector<Particle>::size_type i = 0; i < particles.size(); i++) {
            particles[i].resetNormal();
        }

        for(int x=0; x<widthNum - 1; x++) {
            for(int y=0; y<heightNum - 1; y++) {
                glm::vec3 normal = getNormal(getParticle(x+1,y),getParticle(x,y),getParticle(x,y+1));
                getParticle(x+1,y)->addNormal(normal);
                getParticle(x,y)->addNormal(normal);
                getParticle(x,y+1)->addNormal(normal);

                normal = getNormal(getParticle(x+1,y+1),getParticle(x+1,y),getParticle(x,y+1));
                getParticle(x+1,y+1)->addNormal(normal);
                getParticle(x+1,y)->addNormal(normal);
                getParticle(x,y+1)->addNormal(normal);

                normal = getNormal(getParticle(x+1, y+1),getParticle(x+1, y),getParticle(x,y));
                getParticle(x+1,y+1)->addNormal(normal);
                getParticle(x+1,y)->addNormal(normal);
                getParticle(x,y)->addNormal(normal);

                normal = getNormal(getParticle(x+1, y+1),getParticle(x,y),getParticle(x, y+1));
                getParticle(x+1,y+1)->addNormal(normal);
                getParticle(x,y)->addNormal(normal);
                getParticle(x,y+1)->addNormal(normal);
            }
        }
        int p;
        glBegin(GL_TRIANGLES);
        for(int x=0; x<widthNum - 1; x++)
        {
            p = x % 2;
            for(int y=0; y<heightNum - 1; y++)
            {
                glm::vec3 color;
                if ((y+p)%2) {
                    color = glm::vec3(0.082f,0.48f,0.988f);
                }
                else {
                    color = glm::vec3(0.96f,0.75f,0.572f);
                }

                if ((x+y)%2 == 0){
                    drawTriangle(getParticle(x+1,y), getParticle(x,y), getParticle(x,y+1), color);
                    drawTriangle(getParticle(x+1,y+1), getParticle(x+1,y), getParticle(x,y+1), color);
                } else {
                    drawTriangle(getParticle(x+1,y+1), getParticle(x+1,y), getParticle(x,y), color);
                    drawTriangle(getParticle(x+1,y+1), getParticle(x,y), getParticle(x,y+1), color);
                }
            }
        }
        glEnd();
    }

    void drawBalls(bool showForce, float size) {
        glMatrixMode(GL_MODELVIEW);
        glColor3f(1.0f, 1.0f, 1.0f);
        for(std::vector<Particle>::size_type i = 0; i < particles.size(); i++) {
            glPushMatrix();
            glTranslatef(particles[i].pos.x, particles[i].pos.y, particles[i].pos.z);
            if (showForce) {glColor3f(particles[i].color.x, particles[i].color.y, particles[i].color.z);}
            glutSolidSphere(size, 5, 5);
            glPopMatrix();
        }
    }

    // For .obj file output.
    void writeOBJ(ofstream& output) {

        #pragma omp parallel for
        for(std::vector<Particle>::size_type i = 0; i < particles.size(); i++) {
            particles[i].resetNormal();
        }

        for(int x=0; x<widthNum - 1; x++) {
            for(int y=0; y<heightNum - 1; y++) {
                glm::vec3 normal = getNormal(getParticle(x+1,y),getParticle(x,y),getParticle(x,y+1));
                getParticle(x+1,y)->addNormal(normal);
                getParticle(x,y)->addNormal(normal);
                getParticle(x,y+1)->addNormal(normal);

                normal = getNormal(getParticle(x+1,y+1),getParticle(x+1,y),getParticle(x,y+1));
                getParticle(x+1,y+1)->addNormal(normal);
                getParticle(x+1,y)->addNormal(normal);
                getParticle(x,y+1)->addNormal(normal);

                normal = getNormal(getParticle(x+1, y+1),getParticle(x+1, y),getParticle(x,y));
                getParticle(x+1,y+1)->addNormal(normal);
                getParticle(x+1,y)->addNormal(normal);
                getParticle(x,y)->addNormal(normal);

                normal = getNormal(getParticle(x+1, y+1),getParticle(x,y),getParticle(x, y+1));
                getParticle(x+1,y+1)->addNormal(normal);
                getParticle(x,y)->addNormal(normal);
                getParticle(x,y+1)->addNormal(normal);
            }
        }

        output << "nv " << numOfParticles << endl;
        output << "nvn " << numOfParticles << endl;
        
        for(std::vector<Particle>::size_type i = 0; i < particles.size(); i++) {
            output << "v " << particles[i].pos.x << " " << particles[i].pos.y << " " << particles[i].pos.z << endl;
        }
        
        for(std::vector<Particle>::size_type i = 0; i < particles.size(); i++) {
            output << "vn " << particles[i].normal.x << " " << particles[i].normal.y << " " << particles[i].normal.z << endl;
        }

        int i1, i2, i3;
        for(int x=0; x<widthNum - 1; x++)
        {
            for(int y=0; y<heightNum - 1; y++)
            {
                if ((x+y)%2 == 0){
                    i1 = getParticle(x+1,y)->index; 
                    i2 = getParticle(x,y)->index;
                    i3 = getParticle(x,y+1)->index;
                    output << "f " << i1 << "//" << i1 << " " << i2 << "//" << i2 << " " << i3 << "//" << i3 << endl;
                    
                    i1 = getParticle(x+1,y+1)->index; 
                    i2 = getParticle(x+1,y)->index; 
                    i3 = getParticle(x,y+1)->index;
                    output << "f " << i1 << "//" << i1 << " " << i2 << "//" << i2 << " " << i3 << "//" << i3 << endl;

                } else {
                    i1 = getParticle(x+1,y+1)->index;
                    i2 = getParticle(x+1,y)->index;
                    i3 = getParticle(x,y)->index;
                    output << "f " << i1 << "//" << i1 << " " << i2 << "//" << i2 << " " << i3 << "//" << i3 << endl;

                    i1 = getParticle(x+1,y+1)->index;
                    i2 = getParticle(x,y)->index;
                    i3 = getParticle(x,y+1)->index;
                    output << "f " << i1 << "//" << i1 << " " << i2 << "//" << i2 << " " << i3 << "//" << i3 << endl;
                }
            }
        }        

    }

    // For collisions.
    void withBall(glm::vec3 c, float r, glm::vec3 s) {
        ballCol = true;
        center = c; radius = r;
        speed = s;
    }

    void withCube(glm::vec3 in, glm::vec3 ax) {
        cubeCol = true;
        cubeMin = in; cubeMax = ax;
    }

    // For timeStep.
    void timeStep() {

        //Aerodynamic force.
        windForce();

        // Triangle mesh force.
        for(std::vector<TriangleMesh>::size_type i = 0; i < tms.size(); i++) {
            tms[i].timeStep();
        }

        //Edge Spring force.
        for(std::vector<EdgeSpring>::size_type i = 0; i < ess.size(); i++) {
            ess[i].timeStep();
        }

        // Calculate cube info for late use.
        glm::vec3 v1 = glm::vec3(cubeMin.x, cubeMin.y, cubeMin.z);
        glm::vec3 v2 = glm::vec3(cubeMin.x, cubeMin.y, cubeMax.z);
        glm::vec3 v3 = glm::vec3(cubeMax.x, cubeMin.y, cubeMax.z);
        glm::vec3 v4 = glm::vec3(cubeMax.x, cubeMin.y, cubeMin.z);

        glm::vec3 v5 = glm::vec3(cubeMin.x, cubeMax.y, cubeMax.z);
        glm::vec3 v6 = glm::vec3(cubeMin.x, cubeMax.y, cubeMin.z);
        glm::vec3 v7 = glm::vec3(cubeMax.x, cubeMax.y, cubeMin.z);
        glm::vec3 v8 = glm::vec3(cubeMax.x, cubeMax.y, cubeMax.z);

        //Update particles positions with collision detection.
        #pragma omp parallel for
        for(std::vector<Particle>::size_type i = 0; i < particles.size(); i++) {
            glm::vec3 old = particles[i].pos;
            particles[i].timeStep();
            float intersect;

            // Adjust for floor collision.
            if (hasFloor) {
                if (hitFloor(old, particles[i].pos, intersect)) {
                    particles[i].pos = old + (particles[i].pos-old)*intersect;
                    // Set velocity to 0 to make energy conserved.
                    particles[i].vel = glm::vec3(0.0f, 0.0f, 0.0f);
                }
            }

            // Adjust for ball collision.
            if (ballCol) {
                glm::vec3 diff = particles[i].pos - center;
                float diffLength = glm::length(diff);
                if (diffLength < radius) {

                    //Naive way, move particle to the surface in the direction of the normal.
                    particles[i].move(glm::normalize(diff)*(radius-diffLength));
                    
                    //Advanced way (still testing, stickness weird), calculate the intersection and move particle to the intersection.
                    //glm::vec3 backed = particles[i].pos - speed;
                    //ballIntersect(old, backed, intersect); 
                    //particles[i].pos = old + (backed-old)*intersect;
                    //particles[i].move(speed);

                    // Set velocity to 0 to make energy conserved.
                    particles[i].vel = glm::vec3(0.0f, 0.0f, 0.0f);
                }
            }

            // Adjust for cube collision.
            if (cubeCol) {
                bool insideX = particles[i].pos.x < cubeMax.x and particles[i].pos.x > cubeMin.x;
                bool insideY = particles[i].pos.y < cubeMax.y and particles[i].pos.y > cubeMin.y;
                bool insideZ = particles[i].pos.z < cubeMax.z and particles[i].pos.z > cubeMin.z;

                if (insideX and insideY and insideZ) {
                    if (hitSurface(old, particles[i].pos, intersect, v1, v2, v5, v6)) { //Left.
                        particles[i].pos = old + (particles[i].pos-old)*intersect;
                    } else if (hitSurface(old, particles[i].pos, intersect, v2, v3, v8, v5)) { //Front.
                        particles[i].pos = old + (particles[i].pos-old)*intersect;
                    } else if (hitSurface(old, particles[i].pos, intersect, v3, v4, v7, v8)) { //Right.
                        particles[i].pos = old + (particles[i].pos-old)*intersect;
                    } else if (hitSurface(old, particles[i].pos, intersect, v1, v6, v7, v4)) { //Back.
                        particles[i].pos = old + (particles[i].pos-old)*intersect;
                    } else if (hitSurface(old, particles[i].pos, intersect, v5, v8, v7, v6)) { //Top.
                        particles[i].pos = old + (particles[i].pos-old)*intersect;
                    } else if (hitSurface(old, particles[i].pos, intersect, v1, v4, v3, v2)) { //Bottom.
                        particles[i].pos = old + (particles[i].pos-old)*intersect;
                    } else {
                        cout << "This line should never be printed!!!" << endl;
                    }
                    // Set velocity to 0 to make energy conserved.
                    particles[i].vel = glm::vec3(0.0f, 0.0f, 0.0f);
                }
            }
        }
        // Reset in case of scene changes.
        ballCol = false;
        cubeCol = false;
        hasFloor = false;
    }

    // For reset to initial state.
    void reset() {
        particles = particlesBACKUP;
    }
};
