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

Couleur spec_apport(Vector3 pts_espace, Vector3 normal, Lumiere lum, Sphere obj, Vector3 ray, Scene scene, Rayon * reflected){
    Couleur c;
    Vector3 S = reflected_ray(normal, ray);
    Vector3 light_dir;
    light_dir = light_direction(pts_espace, lum.pos);
    light_dir.x = light_dir.x;
    light_dir.y = light_dir.y;
    light_dir.z = light_dir.z;
    c.r = (int)(obj.para.ks * pow(max(dot(S,light_dir),0),0.5) * lum.intensity);
    c.g = (int)(obj.para.ks * pow(max(dot(S,light_dir),0),0.5)* lum.intensity);
    c.b = (int)(obj.para.ks * pow(max(dot(S,light_dir),0),0.5) * lum.intensity);

    Rayon ray1;
    ray1.u = S;
    ray1.origine = pts_espace;
    *reflected = ray1;
    return c;
}