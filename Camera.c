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