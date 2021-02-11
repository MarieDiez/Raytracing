
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

int N=4;
int N_L=2;
int level = 5;

typedef struct _Couleur {
  int r, g, b;
} Couleur;

typedef struct _Vector3{
	float x;
	float y;
    float z;
} Vector3;

typedef struct _Lumiere{
    Vector3 pos;
    int intensity;
    Couleur light_color;
} Lumiere;

typedef struct _Image{
	int width;
	int height;
    Couleur ** pixels;
} Image;

typedef struct _Material_para{
    float kd_intensity;
    Couleur kd;
    float ks;
} Material_para;

typedef struct _Sphere{
    Material_para para;
    Vector3 pos;
    char * name;
    Couleur c;
    int r;
} Sphere;

typedef struct _Camera{
    float f;
    Vector3 origine;
    Vector3 up;
    Vector3 right;
    Vector3 front;
} Camera;

typedef struct _Scene{
	Sphere * list_obj;
    Lumiere * list_lum;
    Camera camera; 
} Scene;

typedef struct _Rayon{
    Vector3 u;
    Vector3  origine;
} Rayon;

typedef struct _Light{
    Rayon light_ray;
    Couleur  light_color;
} Light;


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

int save_image_as_ppm(Image I, const char* fichier){
	FILE* F = fopen(fichier,"w+");
	if (!F){
		return -1;
    }
	fprintf(F,"P3\n%d %d\n255\n",I.width,I.height);
	for(int j=0;j<I.height;j++){
        for(int i=0;i<I.width;i++){
		    fprintf(F,"%d %d %d ",I.pixels[i][j].r,I.pixels[i][j].g,I.pixels[i][j].b);
        }
    }
	fclose(F);
	return 0;
}

float dot(Vector3 v, Vector3 u) {
    float result = 0.0;
    result += v.x*u.x;
    result += v.y*u.y;
    result += v.z*u.z;
    return result;
}

Rayon pixelToRay(float pixelX, float pixelY, Camera camera){
    Rayon rayon;
    rayon.origine = camera.origine;
    Vector3 pixelCenter;
    pixelCenter.x = pixelX;
    pixelCenter.y = pixelY;
    pixelCenter.z = -1;
    rayon.u = normalize(pixelCenter);
    return rayon;
}

Vector3 light_direction(Vector3 pts_espace, Vector3 lum_pos){
    Vector3 direction;
    direction.x = -pts_espace.x + lum_pos.x;
    direction.y = -pts_espace.y + lum_pos.y;
    direction.z = -pts_espace.z +  lum_pos.z;
    direction = normalize(direction);
    return direction;
}

float max(float x, float y){
    if (x > y)
        return x;
    return y;
}

Couleur diffuse_apport(Vector3 pts_espace, Vector3 normal, Lumiere lum, Sphere obj){
    Couleur c;
    Vector3 lumv1;
    lumv1.x = lum.pos.x;
    lumv1.y = lum.pos.y;
    lumv1.z = lum.pos.z;
    c.r = (int)(obj.para.kd_intensity * (obj.para.kd.r*lum.light_color.r*255) * max(dot(normal,light_direction(pts_espace, lumv1)),0) * lum.intensity);
    c.g = (int)(obj.para.kd_intensity * (obj.para.kd.g*lum.light_color.g*255) * max(dot(normal,light_direction(pts_espace, lumv1)),0) * lum.intensity);
    c.b = (int)(obj.para.kd_intensity* (obj.para.kd.b*lum.light_color.b*255) * max(dot(normal,light_direction(pts_espace, lumv1)),0) * lum.intensity);
    return c;
}

Vector3 reflected_ray(Vector3 normal, Vector3 incident){
    Vector3 reflected;
    float c1 = dot( normal, incident );
    reflected.x = incident.x - (2 * normal.x * c1 );
    reflected.y = incident.y - (2 * normal.y * c1 );
    reflected.z = incident.z - (2 * normal.z * c1 );
    return normalize(reflected);
}

bool intersect_obj(Sphere obj, Rayon ray, float * t1, float * t2, Vector3 * normal, Vector3 * P){
    Vector3 o_minus_c;
    o_minus_c.x = ray.origine.x - obj.pos.x;
    o_minus_c.y = ray.origine.y - obj.pos.y;
    o_minus_c.z = ray.origine.z - obj.pos.z;
    float b = 2*dot(ray.u,o_minus_c); 
    float c = dot(o_minus_c,o_minus_c) - (obj.r*obj.r);
    float delta = b*b - 4*c;
    if (delta < 0)
        return 0;
    *t1 = (-b - sqrt(delta)) / 2;
    *t2 = (-b + sqrt(delta)) / 2;
    Vector3 p1;
    p1.x = ray.u.x * (*t1) + ray.origine.x;
    p1.y = ray.u.y * (*t1) + ray.origine.y;
    p1.z = ray.u.z * (*t1) + ray.origine.z;
    // meme direction
    Vector3 dir =light_direction(p1, ray.origine);
    if (dot(dir, ray.u) > 0){
        return 0;
    }
    Vector3 p_minus_c;
    p_minus_c.x = p1.x -obj.pos.x;
    p_minus_c.y = p1.y - obj.pos.y;
    p_minus_c.z = p1.z - obj.pos.z;
    *normal = normalize(p_minus_c);
    *P = p1;
    return 1;
}

bool intersect_obj_bool(Sphere obj, Rayon ray, int l){
   Vector3 o_minus_c;
    o_minus_c.x = ray.origine.x - obj.pos.x;
    o_minus_c.y = ray.origine.y - obj.pos.y;
    o_minus_c.z = ray.origine.z - obj.pos.z;
    float b = 2*dot(ray.u,o_minus_c); 
    float c = dot(o_minus_c,o_minus_c) - (obj.r*obj.r);
    float delta = b*b - 4*c;
    if (delta < 0)
        return 0;
    return 1;
}

int visible(Lumiere lum, Vector3 pts_espace, Scene scene, int l){
    Vector3 lum_pts_dir = light_direction(pts_espace, lum.pos);
    Rayon ray;
    ray.u = lum_pts_dir,
    ray.origine = lum.pos;
    for(int k=0;k<N;k++){
            Sphere obj = scene.list_obj[k];
            if (intersect_obj_bool(obj, ray, l)){
                if (l == k){
                    return 1;
                } else {
                    return 0;
                }
            }
    }
    return 0;
}

Couleur spec_apport(Vector3 pts_espace, Vector3 normal, Lumiere lum, Sphere obj, Vector3 ray, Scene scene, Rayon * reflected){
    Couleur c;
    Vector3 S = reflected_ray(normal, ray);
    Vector3 light_dir;
    light_dir = light_direction(pts_espace, lum.pos);
    light_dir.x = light_dir.x;
    light_dir.y = light_dir.y;
    light_dir.z = light_dir.z;
    c.r = (int)(obj.para.ks * pow(max(dot(S,light_dir),0),5) * lum.intensity);
    c.g = (int)(obj.para.ks * pow(max(dot(S,light_dir),0),5)* lum.intensity);
    c.b = (int)(obj.para.ks * pow(max(dot(S,light_dir),0),5) * lum.intensity);

    Rayon ray1;
    ray1.u = S;
    ray1.origine = pts_espace;
    *reflected = ray1;
    return c;
}

float deg2rad(float deg) 
{ return deg * M_PI / 180; } 

Couleur lancer_rayon2(Rayon ray, Scene scene, Image img, int colorR,int colorG,int colorB, int level_rec, int indexK){
    Couleur rec;
    rec.r= 0;
    rec.g= 0;
    rec.b= 0;
    if (level_rec > 0){
        for(int k=0;k<N;k++){
            if (k != indexK){
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
                            Rayon reflected;
                            Couleur res = diffuse_apport(pts_espace, normal, scene.list_lum[l], obj);
                            Couleur res2 = spec_apport(pts_espace, normal, scene.list_lum[l], obj, ray.u, scene, &reflected);
                            rec = lancer_rayon2(reflected, scene, img, colorR,colorG,colorB,level_rec-1,k);
                            if( res.r< 0){
                                res.r = 0;
                            }
                            if(res.g  < 0){
                                res.g = 0;
                            }
                            if(res.b< 0){
                                res.b = 0;
                            }
                            if(res2.r< 0){
                                res2.r = 0;
                            }
                            if(res2.g  < 0){
                                res2.g = 0;
                            }
                            if(res2.b< 0){
                                res2.b = 0;
                            }
                            if(rec.r< 0){
                                rec.r = 0;
                            }
                            if(rec.g  < 0){
                                rec.g = 0;
                            }
                            if(rec.b< 0){
                                rec.b = 0;
                            }
                            diffuse.r += res.r + res2.r + (int)(obj.para.ks*rec.r);
                            diffuse.g += res.g+ res2.g+ (int)(obj.para.ks*rec.g);
                            diffuse.b += res.b + res2.b+ (int)(obj.para.ks*rec.b);
                        } else{
                            diffuse.r=0;
                            diffuse.g=0;
                            diffuse.b=0;
                        }
                        colorR =  diffuse.r;
                        colorG =  diffuse.g;
                        colorB = diffuse.b;
                    } 
                } 
            }
        }
        rec.r += colorR;
        rec.g += colorG;
        rec.b += colorB;
    }
    return rec;
}

void lancer_rayon(Rayon ray, Scene scene, Image img,int i, int j, int colorR,int colorG,int colorB){
    Couleur background;
    Couleur rec;
    rec.r = 0;
    rec.g = 0;
    rec.b = 0;
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
                    Rayon reflected;
                    Couleur res = diffuse_apport(pts_espace, normal, scene.list_lum[l], obj);
                    Couleur res2 = spec_apport(pts_espace, normal, scene.list_lum[l], obj, ray.u, scene, &reflected);
                    if  (level > 0){
                        rec = lancer_rayon2(reflected, scene, img,colorR,colorG,colorB,level,k);
                    }
                    if( res.r< 0){
                        res.r = 0;
                    }
                    if(res.g  < 0){
                        res.g = 0;
                    }
                    if(res.b< 0){
                        res.b = 0;
                    }
                    if(res2.r< 0){
                        res2.r = 0;
                    }
                    if(res2.g  < 0){
                        res2.g = 0;
                    }
                    if(res2.b< 0){
                        res2.b =0;
                    }
                    if(rec.r< 0){
                        rec.r = 0;
                    }
                    if(rec.g  < 0){
                        rec.g = 0;
                    }
                    if(rec.b< 0){
                        rec.b = 0;
                    }
                    diffuse.r += res.r + res2.r + (int)(obj.para.ks*rec.r);
                    diffuse.g += res.g+ res2.b+ (int)(obj.para.ks*rec.g);
                    diffuse.b += res.b + res2.b+ (int)(obj.para.ks*rec.b);
                } else{
                    diffuse.r=0;
                    diffuse.g=0;
                    diffuse.b=0;
                }
            } 
            img.pixels[i][j] = diffuse;
        } else {
            if (img.pixels[i][j].r == 0 && img.pixels[i][j].g == 0 && img.pixels[i][j].b == 0){
                background.r = (ray.u.y+1)*255;//255
                background.g = (ray.u.y+1)*255;//255
                background.b = (ray.u.y+1)*0;
                img.pixels[i][j] = background;
            }
        }
    }
}

void imageRender(Scene scene){
    Image img;
    img.width = 500;
    img.height = 500;

    img.pixels = malloc(img.height * (img.width * sizeof(Couleur)));
    for (int i=0; i < img.width; i++){
            img.pixels[i] = malloc(img.height * sizeof(Couleur));
    }
    Couleur sphereColor;
    sphereColor.r = 255;
    sphereColor.g = 255;
    sphereColor.b = 255;

    for(int j=0; j < img.height; j++){
        for(int i=0; i < img.width; i++){
            float fov = 90;
            float scale =  tan(deg2rad(fov* 0.5));
            float Px = (2 * ((i + 0.5) / img.width) - 1) * scale; 
            float Py = (1 - 2 * ((j+ 0.5) / img.height)) * scale; 
            Rayon ray = pixelToRay(Px, Py, scene.camera);
            int colorR=0;
            int colorG=0;
            int colorB=0;
            lancer_rayon(ray, scene,img,i,j,colorR,colorG,colorB);
        }
    }
    save_image_as_ppm(img, "image.ppm");
}

int main(){
    Camera  camera;
    Vector3 pos_camera;
    pos_camera.x = 0; // <-
    pos_camera.y = 0; // bas
    pos_camera.z = 0;
    camera.origine = pos_camera;
    camera.f = -1;

    Sphere sphere1;
    Couleur c1;
    c1.r = 1;
    c1.g = 0;
    c1.b = 0;
    Vector3 sphere_pos;
    sphere_pos.x = 5.6;
    sphere_pos.y = 1;
    sphere_pos.z = -10;
    Material_para para;
    para.ks = 1;
    para.kd_intensity = 0.7;
    para.kd = c1;
    sphere1.para = para;
    sphere1.c = c1;
    sphere1.pos = sphere_pos;
    sphere1.r = 3;
    sphere1.name = "Sphere1";

    Sphere sphere2;
    Couleur c2;
    c2.r = 0;
    c2.g = 1;
    c2.b = 0;
    Vector3 sphere_pos2;
    sphere_pos2.x = -5.6;
    sphere_pos2.y = 1;
    sphere_pos2.z = -10;
    Material_para para2;
    para2.ks = 1;
    para2.kd_intensity = 0.5;
    para2.kd = c2;
    sphere2.para = para2;
    sphere2.c = c2;
    sphere2.pos = sphere_pos2;
    sphere2.r = 3;
    sphere2.name = "Sphere2";

    Sphere sphere3;
    Couleur c3;
    c3.r = 1;
    c3.g = 1;
    c3.b = 1;
    Vector3 sphere_pos3;
    sphere_pos3.x = 0;
    sphere_pos3.y =  0;
    sphere_pos3.z = -8;
    Material_para para3;
    para3.ks = 0.1;
    para3.kd_intensity = 0.05;
    para3.kd = c3;
    sphere3.para = para3;
    sphere3.c = c3;
    sphere3.pos = sphere_pos3;
    sphere3.r = 2;
    sphere3.name = "Sphere3";

    Sphere sol;
    Couleur c4;
    c4.r = 1;
    c4.g = 0;
    c4.b = 1;
    Vector3 sol_pos;
    sol_pos.x = 0;
    sol_pos.y =  -22;
    sol_pos.z = -12;
    Material_para para4;
    para4.ks = 0.1;
    para4.kd_intensity = 0.05;
    para4.kd = c4;
    sol.para = para4;
    sol.pos = sol_pos;
    sol.r = 16;
    sol.name = "sol";

    Lumiere lum;
    Vector3 position_lum;
    position_lum.x = 0;
    position_lum.y = 0;
    position_lum.z = 0;
    lum.pos = position_lum;
    lum.intensity = 10;
    Couleur lightColor;
    lightColor.r = 1;
    lightColor.g = 1;
    lightColor.b = 1;
    lum.light_color = lightColor;

    Lumiere lum2;
    Vector3 position_lum2;
    position_lum2.x = 0;
    position_lum2.y = 10;
    position_lum2.z = -15;
    lum2.pos = position_lum2;
    lum2.intensity = 10;
    Couleur lightColor2;
    lightColor2.r = 1;
    lightColor2.g = 1;
    lightColor2.b = 1;
    lum2.light_color = lightColor2;

    Scene scene;
    scene.camera = camera;
    scene.list_obj = malloc(N* sizeof(Sphere));
    scene.list_obj[0] = sphere2;
    scene.list_obj[1] = sphere1;
    scene.list_obj[2] = sphere3;
    scene.list_obj[3] = sol;

    scene.list_lum = malloc(N_L* sizeof(Lumiere));
    scene.list_lum[0] = lum;
    scene.list_lum[1] = lum2;
    
    imageRender(scene);
   
    return EXIT_SUCCESS;
}