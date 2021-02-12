float deg2rad(float deg) 
{ return deg * M_PI / 180; } 

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