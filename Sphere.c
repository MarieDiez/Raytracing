int N=4;

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