#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "Couleur.h"
#include "Vector3.h"
#include "Lumiere.h"
#include "Image.h"
#include "Material_para.h"
#include "Sphere.h"
#include "Camera.h"
#include "Scene.h"
#include "Rayon.h"

#include "Vector3.c"
#include "Image.c"
#include "Camera.c"
#include "Rayon.c"
#include "Couleur.c"
#include "Sphere.c"

int N_L=3;
int level = 5;

/* A FAIRE EN PLUS
bool hitGround(Rayon r, float floorPos){
  float dist = -(r.origine.y + floorPos)/r.u.y;
  if(dist > 0.0001 && dist < 100)
    return true;
  return false;
}*/

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
                }
            } 
            img.pixels[i][j] = diffuse;
        } else {
            if (img.pixels[i][j].r == 0 && img.pixels[i][j].g == 0 && img.pixels[i][j].b == 0){
                background.r = (ray.u.y+1)*0;//255
                background.g = (ray.u.y+1)*0;//255
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
    pos_camera.x = 0;
    pos_camera.y = 0; 
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
    para.ks = 0.2;
    para.kd_intensity = 0.15;
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
    para2.ks = 0.2;
    para2.kd_intensity = 0.15;
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
    para3.ks = 0.2;
    para3.kd_intensity = 0.12;
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
    para4.ks = 0.2;
    para4.kd_intensity = 0.05;
    para4.kd = c4;
    sol.para = para4;
    sol.pos = sol_pos;
    sol.r = 16;
    sol.name = "sol";

    Lumiere lum;
    Vector3 position_lum;
    position_lum.x = 0;
    position_lum.y = 22;
    position_lum.z = -20;
    lum.pos = position_lum;
    lum.intensity = 6;
    Couleur lightColor;
    lightColor.r = 1;
    lightColor.g = 1;
    lightColor.b = 1;
    lum.light_color = lightColor;


    Lumiere lum1;
    Vector3 position_lum1;
    position_lum1.x = 0;
    position_lum1.y = 22;
    position_lum1.z = -20;
    lum1.pos = position_lum1;
    lum1.intensity = 2;
    Couleur lightColor1;
    lightColor1.r = 2;
    lightColor1.g = 1;
    lightColor1.b = 0;
    lum1.light_color = lightColor1;

    Lumiere lum2;
    Vector3 position_lum2;
    position_lum2.x = 0;
    position_lum2.y = 0;
    position_lum2.z = 0;
    lum2.pos = position_lum2;
    lum2.intensity = 5;
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
    scene.list_lum[2] = lum1;
    
    imageRender(scene);
   
    return EXIT_SUCCESS;
}