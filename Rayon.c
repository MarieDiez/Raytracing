Vector3 light_direction(Vector3 pts_espace, Vector3 lum_pos){
    Vector3 direction;
    direction.x = -pts_espace.x + lum_pos.x;
    direction.y = -pts_espace.y + lum_pos.y;
    direction.z = -pts_espace.z +  lum_pos.z;
    direction = normalize(direction);
    return direction;
}
Vector3 reflected_ray(Vector3 normal, Vector3 incident){
    Vector3 reflected;
    float c1 = dot( normal, incident );
    reflected.x = incident.x - (2 * normal.x * c1 );
    reflected.y = incident.y - (2 * normal.y * c1 );
    reflected.z = incident.z - (2 * normal.z * c1 );
    return normalize(reflected);
}