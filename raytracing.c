#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// faire fichier avec leurs fonctions

int N=2;
int N_L=2;

/* STRUCTURE */
typedef struct _Couleur {
  int r, g, b;   
} Couleur;

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
    Vector3 pos;
    int intensity;
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

typedef struct _Sphere{
    Material_para para;
    Vector3 pos;
    Couleur c;
    int r;
} Sphere;

typedef struct _Camera{
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
    Vector3 u;
    Vector3  point;
} Rayon;

/* FONCTIONS */
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

Rayon pixelToRay(int pixelX, int pixelY, Camera camera){
    /*float scale = tan(90/180.0);
    float x= (2.0 * pixelX - 1.0);
    float y= (2.0 * pixelY - 1.0);
    x *= 1 * scale;
    y *= scale;*/
    Rayon rayon;
    rayon.point = camera.c;
    Vector3 pixelCenter;
    float ordo = 1;
    float abs = 1;
    pixelCenter.x = pixelX - camera.p.x;
    pixelCenter.y = pixelY - camera.p.y;
    pixelCenter.z = - camera.f;
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

Couleur diffuse_apport(Vector3 pts_espace, Vector3 normal, Lumiere lum, Sphere obj){
    Couleur c;
    Vector3 lumv1;
    lumv1.x = -lum.pos.x;
    lumv1.y = -lum.pos.y;
    lumv1.z = -lum.pos.z;

    c.r = obj.para.kd * obj.c.r * dot(normal,light_direction(pts_espace, lumv1)) * lum.intensity;
    c.g = obj.para.kd * obj.c.g * dot(normal,light_direction(pts_espace, lumv1)) * lum.intensity;
    c.b = obj.para.kd * obj.c.b * dot(normal,light_direction(pts_espace, lumv1)) * lum.intensity;
    return c;
}

Vector3 reflected_ray(Vector3 normal, Vector3 incident){
    Vector3 reflected;
    float c1 = dot( normal, incident );
    reflected.x = incident.x + (2 * incident.x * c1 );
    reflected.y = incident.y + (2 * incident.y * c1 );
    reflected.z = incident.z + (2 * incident.z * c1 );
    return reflected;
}

Couleur spec_apport(Vector3 pts_espace, Vector3 normal, Lumiere lum, Sphere obj, Vector3 ray){
    Couleur c;
    Vector3 S = reflected_ray(normal, ray);
    c.r = obj.para.ks * pow(dot(S,light_direction(pts_espace, lum.pos)),3) * lum.intensity;
    c.g = obj.para.kd * pow(dot(S,light_direction(pts_espace, lum.pos)),3)* lum.intensity;
    c.b = obj.para.kd * pow(dot(S,light_direction(pts_espace, lum.pos)),3) * lum.intensity;
    return c;
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

float deg2rad(float deg) { 
    return deg * M_PI / 180; 
} 

void imageRender(Scene scene){
    Image img;
    img.width = 500;
    img.height = 500;
    float fov = deg2rad(90);

    img.pixels = malloc(img.height * (img.width * sizeof(Couleur)));
    for (int i=0; i < img.width; i++){
            img.pixels[i] = malloc(img.height * sizeof(Couleur));
    }
    Couleur noir;
    noir.r = 0;
    noir.g = 0;
    noir.b = 0;
    for(int k=0;k<N;k++){
        Sphere obj = scene.list_obj[k];
        for(int i=0; i < img.width; i++){
            for(int j=0; j < img.height; j++){
                // Convertir les cordonnées de l'image dans l'espace...
                
                //float imageAspectRatio = img.width / (float)img.height; 
                //float Px = (2 * ((i + 0.5) / img.width) - 1) * tan(fov / 2 * M_PI / 180) * imageAspectRatio; 
                //float Py = (1 - 2 * ((j+ 0.5) / img.height)) * tan(fov / 2 * M_PI / 180); 
                // Rayon ray = pixelToRay(Px, Py, scene.camera);

                // centrer image
                float x= i-img.width/2;
                float y= j-img.height/2;
                Rayon ray = pixelToRay(x, y, scene.camera);
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
                            Couleur res2 = spec_apport(pts_espace, normal, scene.list_lum[l], obj, ray.u);
                            diffuse.r += (res.r + res2.r);
                            diffuse.g += (res.g + res2.g);
                            diffuse.b += (res.b + res2.b);
                            img.pixels[i][j] = diffuse;
                        }
                    }
                } else{
                        if (img.pixels[i][j].r == 0 && img.pixels[i][j].g == 0 && img.pixels[i][j].b == 0){
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


    Sphere sphere1;
    Couleur c1;
    c1.r = 255;
    c1.g = 0;
    c1.b = 0;
    Vector3 sphere_pos;
    sphere_pos.x = 0;
    sphere_pos.y = 0;
    sphere_pos.z = 300;
    Material_para para;
    para.ks = 0.2;
    para.kd = 0.2;
    sphere1.para = para;
    sphere1.c = c1;
    sphere1.pos = sphere_pos;
    sphere1.r = 250;

    Sphere sphere2;
    Couleur c2;
    c2.r = 0;
    c2.g = 255;
    c2.b = 0;
    Vector3 sphere_pos2;
    sphere_pos2.x = 500;
    sphere_pos2.y = 200;
    sphere_pos2.z = 300;
    Material_para para2;
    para2.ks = 0.2;
    para2.kd = 0.2;
    sphere2.para = para2;
    sphere2.c = c2;
    sphere2.pos = sphere_pos2;
    sphere2.r = 150;

    scene.camera = camera;
    scene.list_obj = malloc(N* sizeof(Sphere));
    scene.list_obj[0] = sphere1;
    scene.list_obj[1] = sphere2;

    Lumiere lum;
    Vector3 position_lum;
    position_lum.x = 100;
    position_lum.y = 100;
    position_lum.z = 0;
    lum.pos = position_lum;
    lum.intensity = 7;

    Lumiere lum2;
    Vector3 position_lum2;
    position_lum2.x = -100;
    position_lum2.y = -100;
    position_lum2.z = 0;
    lum2.pos = position_lum2;
    lum2.intensity = 1;

    scene.list_lum = malloc(N_L* sizeof(Lumiere));
    scene.list_lum[0] = lum;
    scene.list_lum[1] = lum2;

    imageRender(scene);

    return EXIT_SUCCESS;
}