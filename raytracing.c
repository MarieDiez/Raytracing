#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// faire fichier avec leurs fonctions
/* STRUCTURE */
int N=3;
int N_L=2;
int level = 10;

// color
typedef struct _Couleur {
  int r, g, b;   /* rouge, vert, bleu */
} Couleur;

// Vecteurs
typedef struct _Vector2{
	float x;
	float y;
} Vector2;

typedef struct _Vector3{
	float x;
	float y;
    float z;
} Vector3;

typedef struct _Lumiere{
	// vect directeur
    Vector3 pos;
    int intensity;
    //light amount=light color∗light intensity.
} Lumiere;

// image
typedef struct _Image{
	int width;
	int height;
    Couleur ** pixels;
} Image;

typedef struct _Material_para{
    float kd;
    float ks;
} Material_para;

//scene
typedef struct _Sphere{
    Material_para para;
    Vector3 pos;
    Couleur c;
    int r;
} Sphere;

typedef struct _Camera{
	// vect directeur
    float f;
    Vector3 c;
    Vector3 p;
    Vector3 up;
} Camera;
typedef struct _Scene{
	Sphere * list_obj;
    Lumiere * list_lum;
    Camera camera; 
} Scene;

typedef struct _Rayon{
	// vect directeur
    Vector3 u;
    Vector3  point;
} Rayon;


Couleur spec_apport(Vector3 pts_espace, Vector3 normal, Lumiere lum, Sphere obj, Vector3 ray, Scene scene);

float norme(Vector3 vect){
    float norm;
    float n=0;
    int dim = 3;
    n += vect.x * vect.x;
    n += vect.y * vect.y;
    n += vect.z * vect.z;
    
    norm= sqrt(n);
    return norm;
}

Vector3 normalize(Vector3 vect){
    Vector3 vectRes;
    float norm = norme(vect);
    if(norm != 0){
        vectRes.x = vect.x / norm;
        vectRes.y = vect.y / norm;
        vectRes.z = vect.z / norm;
    }
    return vectRes;
}

Vector3 light_direction(Vector3 pts_espace, Vector3 lum_pos){
    Vector3 direction;
    direction.x = pts_espace.x - lum_pos.x;
    direction.y = pts_espace.y - lum_pos.y;
    direction.z = pts_espace.z - lum_pos.z;
    direction = normalize(direction);
    return direction;
}

int save_image_as_ppm(Image I, const char* fichier){
	FILE* F = fopen(fichier,"w+");
	if (!F){
		return -1;
    }
	fprintf(F,"P3\n%d %d\n255\n",I.width,I.height);
	for(int i=0;i<I.width;i++){
        for(int j=0;j<I.height;j++){
		    fprintf(F,"%d %d %d ",I.pixels[i][j].r,I.pixels[i][j].g,I.pixels[i][j].b);
        }
    }
	fclose(F);
	return 0;
}

Rayon pixelToRay(float pixelX, float pixelY, Camera camera){
    /*float x= (2.0 * pixelX - 1.0);
    float y= (2.0 * pixelY - 1.0);
    x *= 1 * scale;
    y *= scale;*/

    Rayon rayon;
    rayon.point = camera.c;
    Vector3 pixelCenter;
    float right = 1;
    float up = -1;
    pixelCenter.x = pixelX*up;
    pixelCenter.y = pixelY ;
    pixelCenter.z = -1.7;
    /*pixelCenter.x = pixelX;
    pixelCenter.y = pixelY;
    pixelCenter.z = -camera.f;*/
    rayon.u = normalize(pixelCenter);
    return rayon;
}

float dot(Vector3 v, Vector3 u) {
    double result = 0.0;
    result += v.x*u.x;
    result += v.y*u.y;
    result += v.z*u.z;
    return result;
}



bool intersect_obj(Sphere obj, Rayon ray, float * t1, float * t2, Vector3 * normal, Vector3 * P){
    Vector3 o_minus_c;
    Vector3 direct;
    o_minus_c.x = ray.point.x + obj.pos.x;
    o_minus_c.y = ray.point.y + obj.pos.y;
    o_minus_c.z = ray.point.z + obj.pos.z;
    float b = 2*dot(ray.u,o_minus_c); 
    float c = dot(o_minus_c,o_minus_c) - (obj.r*obj.r);
    float delta = b*b - 4*c;
    if (delta < 0){
        return 0;
    } else{
        *t1 = (-b - sqrt(delta)) / 2;
        *t2 = (-b + sqrt(delta)) / 2;
        Vector3 p1;
        float droot = sqrt(delta);
        p1.x = ray.u.x * (*t1) + ray.point.x;
        p1.y = ray.u.y * (*t1) + ray.point.y;
        p1.z = ray.u.z * (*t1) + ray.point.z;
        Vector3 p_minus_c;
        p_minus_c.x = p1.x - obj.pos.x;
        p_minus_c.y = p1.y - obj.pos.y;
        p_minus_c.z = p1.z - obj.pos.z;
        *normal = normalize(p_minus_c);
        *P = p1;
        return 1;
    }
}
bool intersect_obj_bool(Sphere obj, Rayon ray, int l){
    Vector3 o_minus_c;
    o_minus_c.x = ray.point.x + obj.pos.x;
    o_minus_c.y = ray.point.y + obj.pos.y;
    o_minus_c.z = ray.point.z + obj.pos.z;
    float b = 2*dot(ray.u,o_minus_c); 
    float c = dot(o_minus_c,o_minus_c) - (obj.r*obj.r);
    float delta = b*b - 4*c;
    if (delta < 0){
        return 0;
    } else{
        return 1;
    }
}
int visible(Lumiere lum, Vector3 pts_espace, Scene scene, int l){
    Vector3 lum_pts_dir = light_direction(pts_espace, lum.pos);
    Rayon ray;
    ray.u = lum_pts_dir,
    ray.point = lum.pos;
    for(int k=0;k<N;k++){
        if(k!=l){
            Sphere obj = scene.list_obj[k];
            if (intersect_obj_bool(obj, ray, l)){
                return 0;
            }
        }
    }
    return 1;
}
Couleur diffuse_apport(Vector3 pts_espace, Vector3 normal, Lumiere lum, Sphere obj){
 // iD = KD*c*(n.l) * iLI
    Couleur c;
    Vector3 lumv1;
    lumv1.x = lum.pos.x;
    lumv1.y = lum.pos.y;
    lumv1.z = lum.pos.z;

    c.r = obj.para.kd * obj.c.r * dot(normal,light_direction(pts_espace, lumv1)) * lum.intensity;
    c.g = obj.para.kd * obj.c.g * dot(normal,light_direction(pts_espace, lumv1)) * lum.intensity;
    c.b = obj.para.kd * obj.c.b * dot(normal,light_direction(pts_espace, lumv1)) * lum.intensity;
    return c;
}

Vector3 reflected_ray(Vector3 normal, Vector3 incident){
    /*double cosI = -dot(normal, incident);
    double sinT2 = (1 - cosI*cosI);
    double cosT = sqrt(1-sinT2);*/
    Vector3 reflected;
  /*  reflected.x = incident.x + (cosI-cosT) * normal.x;
    reflected.y = incident.y + (cosI-cosT) * normal.y;
    reflected.z = incident.z + (cosI-cosT) * normal.z;
*/
    float c1 = dot( normal, incident );
    reflected.x = incident.x + (2 * incident.x * c1 );
    reflected.y = incident.y + (2 * incident.y * c1 );
    reflected.z = incident.z + (2 * incident.z * c1 );
    return reflected;
}

void reflexion_rec(Rayon ray, float *r, float *g, float *b, Scene scene){
    for(int k=0;k<N;k++){
        Sphere obj = scene.list_obj[k];
        float t1, t2;
        Vector3 normal;
        Vector3 pts_espace;
        if(intersect_obj(obj, ray, &t1, &t2, &normal, &pts_espace)){
            Couleur diffuse;
            diffuse.r=0;
            diffuse.g=0;
            diffuse.b=0;
            for(int l=0;l<N_L;l++){ 
                if (visible(scene.list_lum[l],pts_espace, scene, k)){
                    Couleur res = diffuse_apport(pts_espace, normal, scene.list_lum[l], obj);
                    Couleur res2 = spec_apport(pts_espace, normal, scene.list_lum[l], obj, ray.u,scene);
                    *r += (res.r);
                    *g += (res.g);
                    *b += (res.b);
                }
            }
        } else {
            *r = 0;
            *g = 0;
            *b = 0;
        }
    } 
}

Couleur spec_apport(Vector3 pts_espace, Vector3 normal, Lumiere lum, Sphere obj, Vector3 ray, Scene scene){
    //id = ks * (S.l)^ns*I 
    Couleur c;
    Vector3 S = reflected_ray(normal, ray);
    float r = 0;
    float g = 0;
    float b = 0;
    Rayon ray1;
    ray1.u = S;
    ray1.point = pts_espace;
    if  (level  > 0){
        level-=1;
        reflexion_rec(ray1,&r,&g,&b,scene);
    }
    c.r = obj.para.ks * pow(dot(S,light_direction(pts_espace, lum.pos)),5) * lum.intensity +r;
    c.g = obj.para.ks * pow(dot(S,light_direction(pts_espace, lum.pos)),5)* lum.intensity+g ;
    c.b = obj.para.ks * pow(dot(S,light_direction(pts_espace, lum.pos)),5) * lum.intensity+b;
    return c;
}


float deg2rad(float deg) 
{ return deg * M_PI / 180; } 



void imageRender(Scene scene){
    Image img;
    img.width = 500;
    img.height = 500;

    img.pixels = malloc(img.height * (img.width * sizeof(Couleur)));
    for (int i=0; i < img.width; i++){
            img.pixels[i] = malloc(img.height * sizeof(Couleur));
    }
    Couleur noir;
    noir.r = 255;
    noir.g = 255;
    noir.b = 0;
    for(int i=0; i < img.width; i++){
        for(int j=0; j < img.height; j++){
            // Convertir les cordonnées de l'image dans l'espace...
            float fov = 90;
            float imageAspectRatio = img.width / (float)img.height; 
            float Px = (2 * ((i + 0.5) / img.width) - 1) * tan(fov / 2 * M_PI / 180) * imageAspectRatio; 
            float Py = (1 - 2 * ((j+ 0.5) / img.height)) * tan(fov / 2 * M_PI / 180); 
            Rayon ray = pixelToRay(Px, Py, scene.camera);
            /*float r = 0;
            float g = 0;
            float b = 0;
            reflexion_rec(ray, &r, &g, &b, scene);
            Couleur col;
            col.r=r;
            col.g=g;
            col.b=b;
            img.pixels[i][j] = col;*/
            for(int k=0;k<N;k++){
                Sphere obj = scene.list_obj[k];
                float t1, t2;
                Vector3 normal;
                Vector3 pts_espace;
                if(intersect_obj(obj, ray, &t1, &t2, &normal, &pts_espace)){
                    Couleur diffuse;
                    diffuse.r=0;
                    diffuse.g=0;
                    diffuse.b=0;
                    for(int l=0;l<N_L;l++){ 
                        if (visible(scene.list_lum[l],pts_espace, scene, k)){
                            Couleur res = diffuse_apport(pts_espace, normal, scene.list_lum[l], obj);
                            Couleur res2 = spec_apport(pts_espace, normal, scene.list_lum[l], obj, ray.u, scene);
                            diffuse.r += (res.r+res2.r);
                            diffuse.g += (res.g+res2.g);
                            diffuse.b += (res.b+res2.b);
                            img.pixels[i][j] = diffuse;
                        }
                    }
                } else{
                        if (img.pixels[i][j].r == 0 && img.pixels[i][j].g == 0 && img.pixels[i][j].b == 0){
                            noir.r = (ray.u.y+1)*255;
                            noir.g = (ray.u.y+1)*255;
                            noir.b = (ray.u.y+1)*0;
                            img.pixels[i][j] = noir;
                        }
                    }
            }
        }
    }
    save_image_as_ppm(img, "image.ppm");
}

int main(){
    Scene scene;
    Camera  camera;
    camera.f = 70;
    Vector3 pos_camera;
    pos_camera.x = 0;
    pos_camera.y = 0;
    pos_camera.z = 0;
    camera.c = pos_camera;

    Couleur c1;
    c1.r = 255;
    c1.g = 0;
    c1.b = 0;
    Sphere sphere1;
    Vector3 sphere_pos;
    sphere_pos.x = 0;
    sphere_pos.y = 0;
    sphere_pos.z = 20;
    Material_para para;
    para.ks = 0.1;
    para.kd = 0.2;
    sphere1.para = para;
    sphere1.c = c1;
    sphere1.pos = sphere_pos;
    sphere1.r = 7;

    Couleur c2;
    c2.r = 0;
    c2.g = 255;
    c2.b = 0;

    Sphere sphere2;
    Vector3 sphere_pos2;
    sphere_pos2.x = -5;
    sphere_pos2.y =  -5;
    sphere_pos2.z = 12;
    Material_para para2;
    para2.ks = 0.2;
    para2.kd = 0.06;
    sphere2.para = para2;
    sphere2.c = c2;
    sphere2.pos = sphere_pos2;
    sphere2.r = 3;

    Couleur c3;
    c3.r = 0;
    c3.g = 0;
    c3.b = 255;
    Sphere sol;
    Vector3 sol_pos;
    sol_pos.x = 90;
    sol_pos.y =  0;
    sol_pos.z = 70;
    Material_para para3;
    para3.ks = 0.15;
    para3.kd = 0.06;
    sol.para = para3;
    sol.c = c3;
    sol.pos = sol_pos;
    sol.r = 70;

    scene.camera = camera;
    scene.list_obj = malloc(N* sizeof(Sphere));
    scene.list_obj[0] = sol;
    scene.list_obj[1] = sphere1;
    scene.list_obj[2] = sphere2;

    Lumiere lum;
    Vector3 position_lum;
    position_lum.x = -30;
    position_lum.y = -60;
    position_lum.z = 0;
    lum.pos = position_lum;
    lum.intensity = 20;
    Lumiere lum2;
    Vector3 position_lum2;
    position_lum2.x = -60;
    position_lum2.y = 100;
    position_lum2.z = 0;
    lum2.pos = position_lum2;
    lum2.intensity = 5;
    scene.list_lum = malloc(N_L* sizeof(Lumiere));
    scene.list_lum[0] = lum;
    scene.list_lum[1] = lum2;

    imageRender(scene);
    
   
    return EXIT_SUCCESS;
}