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

float dot(Vector3 v, Vector3 u) {
    float result = 0.0;
    result += v.x*u.x;
    result += v.y*u.y;
    result += v.z*u.z;
    return result;
}