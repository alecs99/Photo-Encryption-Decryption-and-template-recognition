#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef struct pixel {
    unsigned int blue, green, red;
};
typedef struct detectii{
    double corelatie;
    unsigned int row,column,width,height;
};

void xorshift32(unsigned int **r, unsigned int seed, unsigned int width, unsigned int height) {
    unsigned int i, r0;
    *r = (unsigned int *) malloc(2 * width * height * sizeof(unsigned int));
    r0 = seed;
    for (i = 1; i < 2 * width * height; i++) {
        r0 = r0 ^ r0 << 13;
        r0 = r0 ^ r0 >> 17;
        r0 = r0 ^ r0 << 5;
        (*r)[i] = r0;
    }
}

void incarca_liniarizata(unsigned int **h, struct pixel **a, char *nume_poza) {
    FILE *img;
    unsigned int length, width, x;
    unsigned int long padding;
    int i, j;
    unsigned char b, g, r;
    img = fopen(nume_poza, "rb");
    if (img == NULL) {
        printf("Nu am gasit imaginea");
        return;
    }
    fseek(img, 18, SEEK_SET);
    fread(&width, sizeof(int), 1, img);
    fread(&length, sizeof(int), 1, img);
    fseek(img, 0, SEEK_SET);

    (*h) = malloc(54 * sizeof(unsigned int));
    for (i = 0; i < 54; i++) {
        fread(&x, 1, 1, img);
        (*h)[i] = x;
    }
    fseek(img, 54, SEEK_SET);
    (*a) = malloc(width * length * sizeof(struct pixel));
    if (width % 4 != 0)
        padding = 4 - (3 * width) % 4;
    else padding = 0;
    for (i = length - 1; i >= 0; i--) {
        for (j = 0; j < width; j++) {
            fread(&b, 1, 1, img);
            (*a)[j + i * width].blue = b;
            fread(&g, 1, 1, img);
            (*a)[j + i * width].green = g;
            fread(&r, 1, 1, img);
            (*a)[j + i * width].red = r;
        }
        fseek(img, padding, SEEK_CUR);
    }
    fclose(img);
}

void salveaza_liniarizata(unsigned int *h, struct pixel *a, char *nume_poza) {
    FILE *img;
    unsigned int length, width, j, x, padding, p = 0;
    int i;
    unsigned int b, g, r;
    img = fopen(nume_poza, "wb+");
    if(img==NULL)
    {
        printf("Nu se deschide");
        return;
    }
    for (i = 0; i < 54; i++) {
        x = h[i];
        fwrite(&x, 1, 1, img);
    }
    fseek(img, 18, SEEK_SET);
    fread(&width, sizeof(int), 1, img);
    fread(&length, sizeof(int), 1, img);
    fseek(img, 54, SEEK_SET);
    if (width % 4 != 0)
        padding = 4 - (3 * width) % 4;
    else padding = 0;
    for (i = length - 1; i >= 0; i--) {
        for (j = 0; j < width; j++) {
            b = a[j + i * width].blue;
            fwrite(&b, 1, 1, img);
            g = a[j + i * width].green;
            fwrite(&g, 1, 1, img);
            r = a[j + i * width].red;
            fwrite(&r, 1, 1, img);
        }
        fwrite(&p, 1, padding, img);

    }
    fclose(img);
}

void permutare(unsigned int **p, unsigned int r0, unsigned w, unsigned h) {
    int k;
    unsigned int rev, aux;
    unsigned int *xor;
    xorshift32(&xor, r0, w, h);
    *p = (unsigned int *) malloc(w * h * sizeof(unsigned int));
    for (k = 0; k < w * h; k++) {
        (*p)[k] = k;
    }
    for (k = (w * h) - 1; k >= 0; k--) {
        rev = xor[(w * h) - k] % (k + 1);
        aux = (*p)[k];
        (*p)[k] = (*p)[rev];
        (*p)[rev] = aux;
    }

}

void criptare(char *nume_fisier_initial, char *nume_fisier_criptat, char *cheie, struct pixel *l) {
    FILE *in;
    FILE *cr;
    FILE *ch;
    in = fopen(nume_fisier_initial, "rb");
    cr = fopen(nume_fisier_criptat, "wb");
    ch = fopen(cheie, "r");
    if (in == NULL) {
        printf("Nu se deschide");
        return;
    }
    if (cr == NULL) {
        printf("Nu se deschide");
        return;
    }
    if (ch == NULL) {
        printf("Nu se deschide");
        return;
    }
    int i;
    unsigned int *h;
    unsigned int x;
    h = malloc(54 * sizeof(unsigned int));
    for (i = 0; i < 54; i++) {
        fread(&x, 1, 1, in);
        h[i] = x;
    }
    for (i = 0; i < 54; i++) {
        x = h[i];
        fwrite(&x, 1, 1, cr);
    }
    free(h);
    unsigned int b, g, r;
    unsigned int n, r0, sv, l_img, h_img, padding, pad = 0;
    unsigned int *xor, *p;
    struct pixel *lm;
    int k, j = 0;
    fscanf(ch, "%u", &r0);
    fscanf(ch, "%u", &sv);
    fclose(ch);
    fseek(in, 18, SEEK_SET);
    fread(&l_img, sizeof(unsigned int), 1, in);
    fread(&h_img, sizeof(unsigned int), 1, in);
    fseek(in, 54, SEEK_SET);
    n = l_img * h_img;
    xorshift32(&xor, r0, l_img, h_img);
    lm = malloc(n * sizeof(struct pixel));
    if (l_img % 4 != 0)
        padding = 4 - (3 * l_img) % 4;
    else padding = 0;
    permutare(&p, r0, l_img, h_img);
    for (k = 0; k < n; k++) {
        lm[p[k]].blue = l[k].blue;
        lm[p[k]].green = l[k].green;
        lm[p[k]].red = l[k].red;
    }
    struct pixel *c;
    c = malloc(n * sizeof(struct pixel));
    c[0].blue = (sv & 255) ^ lm[0].blue ^ (xor[n] & 255);
    c[0].green = ((sv >> 8) & 255) ^ lm[0].green ^ ((xor[n] >> 8) & 255);
    c[0].red = (((sv >> 8) >> 8) & 255) ^ lm[0].red ^ (((xor[n] >> 8) >> 8) & 255);
    for (k = 1; k < n; k++) {
        c[k].blue = (c[k - 1].blue ^ lm[k].blue) ^ (xor[n + k] & 255);
        c[k].green = (c[k - 1].green ^ lm[k].green) ^ ((xor[n + k] >> 8) & 255);
        c[k].red = (c[k - 1].red ^ lm[k].red) ^ (((xor[n + k] >> 8) >> 8) & 255);
    }
    for (i = h_img - 1; i >= 0; i--) {
        for (j = 0; j < l_img; j++) {
            b = c[j + i * l_img].blue;
            fwrite(&b, 1, 1, cr);
            g = c[j + i * l_img].green;
            fwrite(&g, 1, 1, cr);
            r = c[j + i * l_img].red;
            fwrite(&r, 1, 1, cr);
        }
        fwrite(&pad, 1, padding, cr);
    }
    free(c);
    free(lm);
    free(xor);
    free(p);
    fclose(in);
    fclose(cr);
}

void permutare_inversa(unsigned int **pinv, unsigned int *p, unsigned width, unsigned height) {
    unsigned int n = width * height;
    int k;
    *pinv = malloc(n * sizeof(unsigned int));
    for (k = 0; k < n; k++) {
        (*pinv)[p[k]] = k;
    }


}

void decriptare(char *nume_fisier_decriptat, char *nume_fisier_criptat, char *nume_fisier_cheie) {
    FILE *out;
    FILE *cr;
    FILE *cheie;
    out = fopen(nume_fisier_decriptat, "wb");
    cr = fopen(nume_fisier_criptat, "rb");
    cheie = fopen(nume_fisier_cheie, "r");
    if (out == NULL) {
        printf("Nu se deschide");
        return;
    }
    if (cr == NULL) {
        printf("Nu se deschide");
        return;
    }
    if (cheie == NULL) {
        printf("Nu se deschide");
        return;
    }
    unsigned int *h;
    unsigned int x;
    int i, k, j;
    h = malloc(54 * sizeof(unsigned int));
    for (i = 0; i < 54; i++) {
        fread(&x, 1, 1, cr);
        h[i] = x;
    }
    fseek(out,0,SEEK_SET);
    for (i = 0; i < 54; i++) {
        x = h[i];
        fwrite(&x, 1, 1, out);
    }
    unsigned int b, g, r;
    unsigned int n, r0, sv, l_img, h_img, padding, pad = 0;
    unsigned int *xor, *p, *pinv;
    struct pixel *enc_liniar, *lm;
    fscanf(cheie, "%u", &r0);
    fscanf(cheie, "%u", &sv);
    fclose(cheie);
    fseek(cr, 18, SEEK_SET);
    fread(&l_img, sizeof(unsigned int), 1, cr);
    fread(&h_img, sizeof(unsigned int), 1, cr);
    fseek(cr, 54, SEEK_SET);
    n = l_img * h_img;
    xorshift32(&xor, r0, l_img, h_img);
    enc_liniar = malloc(n * sizeof(struct pixel));
    lm = malloc(n * sizeof(struct pixel));
    if (l_img % 4 != 0)
        padding = 4 - (3 * l_img) % 4;
    else padding = 0;
    permutare(&p, r0, l_img, h_img);
    permutare_inversa(&pinv, p, l_img, h_img);
    for (i = h_img - 1; i >= 0; i--) {
        for (j = 0; j < l_img; j++) {
            fread(&b, 1, 1, cr);
            enc_liniar[j + i * l_img].blue = b;
            fread(&g, 1, 1, cr);
            enc_liniar[j + i * l_img].green = g;
            fread(&r, 1, 1, cr);
            enc_liniar[j + i * l_img].red = r;
        }
        fseek(cr, padding, SEEK_CUR);

    }
    struct pixel *c;
    c = malloc(n * sizeof(struct pixel));
    c[0].blue = (sv & 255) ^ enc_liniar[0].blue ^ (xor[n] & 255);
    c[0].green = ((sv >> 8) & 255) ^ enc_liniar[0].green ^ ((xor[n] >> 8) & 255);
    c[0].red = (((sv >> 8) >> 8) & 255) ^ enc_liniar[0].red ^ (((xor[n] >> 8) >> 8) & 255);
    for (k = 1; k < n; k++) {
        c[k].blue = (enc_liniar[k - 1].blue ^ enc_liniar[k].blue) ^ (xor[n + k] & 255);
        c[k].green = (enc_liniar[k - 1].green ^ enc_liniar[k].green) ^ ((xor[n + k] >> 8) & 255);
        c[k].red = (enc_liniar[k - 1].red ^ enc_liniar[k].red) ^ (((xor[n + k] >> 8) >> 8) & 255);
    }
    for (k = 0; k < n; k++) {
        lm[pinv[k]].blue = c[k].blue;
        lm[pinv[k]].green = c[k].green;
        lm[pinv[k]].red = c[k].red;

    }
    for (i = h_img - 1; i >= 0; i--) {
        for (j = 0; j < l_img; j++) {
            b = lm[j + i * l_img].blue;
            fwrite(&b, 1, 1, out);
            g = lm[j + i * l_img].green;
            fwrite(&g, 1, 1, out);
            r = lm[j + i * l_img].red;
            fwrite(&r, 1, 1, out);
        }
        fwrite(&pad, 1, padding, out);
    }
    free(lm);
    free(c);
    free(enc_liniar);
    fclose(out);
    fclose(cr);
}

void test_chi(char *nume_fisier, struct pixel *a) {
    FILE *in;
    in = fopen(nume_fisier, "rb");
    unsigned int *r, *g, *b;
    unsigned int l_img, h_img, n, i, j, padding;
    float p;
    float x_red, x_green, x_blue, f, w;
    r = malloc(256 * sizeof(unsigned int));
    g = malloc(256 * sizeof(unsigned int));
    b = malloc(256 * sizeof(unsigned int));
    fseek(in, 18, SEEK_SET);
    fread(&l_img, sizeof(unsigned int), 1, in);
    fread(&h_img, sizeof(unsigned int), 1, in);
    fseek(in, 54, SEEK_SET);
    n = l_img * h_img;
    x_red = 0;
    x_green = 0;
    x_blue = 0;
    if (l_img % 4 != 0)
        padding = 4 - (3 * l_img) % 4;
    else padding = 0;
    for (i = 0; i < 256; i++) {
        r[i] = 0;
        g[i] = 0;
        b[i] = 0;
    }
    for (i = 0; i < h_img; i++) {
        for (j = 0; j < l_img; j++) {
            r[a[i + j * l_img].red] = r[a[i + j * l_img].red] + 1;
            g[a[i + j * l_img].green] = g[a[i + j * l_img].green] + 1;
            b[a[i + j * l_img].blue] = b[a[i + j * l_img].blue] + 1;


        }
        fseek(in, padding, SEEK_CUR);
    }
    f = n / 256.0;
    for (i = 0; i < 256; i++) {
        p = (r[i] - f) * (r[i] - f);
        x_red = x_red + (p / f);
        p = (g[i] - f) * (g[i] - f);
        x_green = x_green + (p / f);
        p = (b[i] - f) * (b[i] - f);
        x_blue = x_blue + (p / f);
    }
    printf("R:%.2f \n", x_red);
    printf("G:%.2f \n", x_green);
    printf("B:%.2f \n", x_blue);
    free(r);
    free(g);
    free(b);
    fclose(in);
}

void grayscale(char *nume_fisier_initial,char *nume_fisier_destinatie)
{
    FILE *in,*out;
    in=fopen(nume_fisier_initial,"rb");
    out=fopen(nume_fisier_destinatie,"wb+");
    if(in==NULL)
    {
        printf("Nu exista");
    }
    unsigned int dim_img, l_img, h_img;
    unsigned char *prgb, header[54], aux;
    prgb=malloc(3*sizeof(unsigned char));
    fseek(in,18,SEEK_SET);
    fread(&l_img, sizeof(unsigned int),1,in);
    fread(&h_img, sizeof(unsigned int),1,in);
    fseek(in,0,SEEK_SET);
    unsigned char c;
    while(fread(&c,1,1,in)==1)
    {
        fwrite(&c,1,1,out);
        fflush(out);
    }
    fclose(in);
    int padding;
    if(l_img % 4 != 0)
        padding=4-(3*l_img)%4;
        else padding=0;
    fseek(out, 54, SEEK_SET);
    int i,j;
    for(i=0;i<h_img;i++)
    {
        for(j=0;j <l_img;j++)
        {
            fread(prgb,3,1,out);
            aux = 0.299*prgb[2] + 0.587*prgb[1] + 0.114*prgb[0];
            prgb[0] = prgb[1] = prgb[2] = aux;
            fseek(out, -3, SEEK_CUR);
            fwrite(prgb,3,1,out);
            fflush(out);
        }
        fseek(out,padding,SEEK_CUR);
    }
    fclose(out);
}
double correlation(struct pixel *fi,struct pixel *s,unsigned int h_s,unsigned int l_s)
{
    double n=h_s*l_s,sigma_fi=0,sigma_s=0,fi_barat=0,s_barat=0,corr=0;
    unsigned int i,j;
    for(i=0;i<n;i++) {
        fi_barat = fi_barat + fi[i].red;
        s_barat = s_barat + s[i].red;
    }
    fi_barat=fi_barat/n;
    s_barat=s_barat/n;
    for(i=0;i<n;i++)
    {
        sigma_fi=sigma_fi+(fi[i].red-fi_barat)*(fi[i].red-fi_barat);
        sigma_s=sigma_s+(s[i].red-s_barat)*(s[i].red-s_barat);
    }
    sigma_fi=sigma_fi/(n-1);
    sigma_fi=sqrt(sigma_fi);
    sigma_s=sigma_s/(n-1);
    sigma_s=sqrt(sigma_s);
    for(i=0;i<n;i++)
    {

        corr=corr+(((fi[i].red-fi_barat)*(s[i].red-s_barat))/(sigma_fi*sigma_s));

    }
    corr=corr/n;
    return corr;


}
void border(char *imagine,unsigned int *prgb, int row, unsigned int column, unsigned int height_sab,
            unsigned int width_sab) {
    FILE *img;
    img = fopen(imagine, "rb");
    unsigned int b,g,r;
    unsigned int width_img, heigth_img;
    int i, j;
    unsigned int long padding,p=0;
    unsigned int x;
    unsigned int *header;
    header=malloc(54*sizeof(unsigned int));
    for (i = 0; i < 54; i++) {
        fread(&x, 1, 1, img);
        header[i]=x;
    }
    fseek(img, 18, SEEK_SET);
    fread(&width_img, sizeof(int), 1, img);
    fread(&heigth_img, sizeof(int), 1, img);
    printf("%u %u \n",width_img,heigth_img);
    struct pixel *im;
    im=malloc(width_img*heigth_img*sizeof(struct pixel));
    if(width_img%4!=0)
        padding=4-(3*width_img)%4;
    else padding=0;
    fseek(img,54,SEEK_SET);
    for(i=heigth_img-1;i>=0;i--)
    {
        for(j=0;j<width_img;j++)
        {
            fread(&b, 1, 1, img);
            im[j + i * width_img].blue = b;
            fread(&g, 1, 1, img);
            im[j + i * width_img].green = g;
            fread(&r, 1, 1, img);
            im[j + i * width_img].red = r;

        }
        fseek(img,padding,SEEK_CUR);
    }
    fclose(img);
    i=heigth_img-row-1;
    printf("h:%u w:%u \n",height_sab,width_sab);
    printf("row:%u col:%u \n",row+1,column+1);
    for(j=column;j<=column+width_sab-1;j++)
    {
        im[j+i*width_img].blue=prgb[0];
        im[j+i*width_img].green=prgb[1];
        im[j+i*width_img].red=prgb[2];
    }
    i=heigth_img-row-height_sab;
    for(j=column;j<=column+width_sab-1;j++)
    {
        im[j+i*width_img].blue=prgb[0];
        im[j+i*width_img].green=prgb[1];
        im[j+i*width_img].red=prgb[2];
    }
    j=column;
    for(i=heigth_img-row-height_sab;i<=heigth_img-row-1;i++)
    {
        im[j+i*width_img].blue=prgb[0];
        im[j+i*width_img].green=prgb[1];
        im[j+i*width_img].red=prgb[2];

    }
    j=column+width_sab-1;
    for(i=heigth_img-row-height_sab;i<=heigth_img-row-1;i++)
    {
        im[j+i*width_img].blue=prgb[0];
        im[j+i*width_img].green=prgb[1];
        im[j+i*width_img].red=prgb[2];

    }
    FILE *out;
    out=fopen(imagine,"wb");
    for (i = 0; i < 54; i++) {
        x=header[i];
        fwrite(&x, 1, 1, out);
    }
    fseek(out,54,SEEK_SET);
    for (i = heigth_img - 1; i >= 0; i--) {

        for (j = 0; j < width_img; j++) {
            b = im[j + i * width_img].blue;
            fwrite(&b, 1, 1, out);
            g = im[j + i * width_img].green;
            fwrite(&g, 1, 1, out);
            r = im[j + i * width_img].red;
            fwrite(&r, 1, 1, out);
        }
        fwrite(&p, 1, padding, out);
    }
    fclose(out);


}
struct detectii *match(char *imagine,char *sablon,double prag,unsigned int *n)
{
    FILE* img;
    FILE *sab;
    img=fopen(imagine,"rb");
    sab=fopen(sablon,"rb");
    if((img==NULL)||(sab==NULL))
    {
        printf("Nu se deschide");
    }
    unsigned int j,y,l_img,h_img,l_s,h_s,m,poz_l,poz_c,w=0;
    unsigned char b,g,r;
    int i,x;
    unsigned int long jump_s,jump_fi,padding_fi,padding_s;
    double corr;
    fseek(img, 18, SEEK_SET);
    fread(&l_img,sizeof(int),1,img);
    fread(&h_img, sizeof(int),1,img);
    fseek(sab, 18, SEEK_SET);
    fread(&l_s,sizeof(int),1,sab);
    fread(&h_s, sizeof(int),1,sab);
    fseek(sab,54,SEEK_SET);
    m=l_s*h_s;
    struct pixel *fi,*s;
    struct detectii *d;
    d=malloc(m*sizeof(struct detectii));
    fi=malloc(m*sizeof(struct pixel));
    s=malloc(m*sizeof(struct pixel));
    if(l_img%4!=0)
        padding_fi=4-(3*l_img)%4;
    else padding_fi=0;
    if(l_s%4!=0)
        padding_s=4-(3*l_s)%4;
    else padding_s=0;
    for(i=h_s-1;i>=0;i--) {

        for (j = 0; j < l_s; j++) {
            fread(&b, 1, 1, sab);
            s[j + i * l_s].blue = b;
            fread(&g, 1, 1, sab);
            s[j + i * l_s].green = g;
            fread(&r, 1, 1, sab);
            s[j + i * l_s].red = r;

        }
        fseek(sab,padding_s,SEEK_CUR);
    }
    printf("Inaltime:%u \n Latime:%u \n",h_s,l_s);
    fseek(img,54,SEEK_SET);
    for (i = 0; i <= h_img-h_s; i++) {
        for (j = 0; j <= l_img-l_s; j++) {
            jump_fi=54+3*(i*l_img+j)+i*padding_fi;
            fseek(img,jump_fi,SEEK_SET);
            for (x = h_s-1; x>=0; x--) {
                for (y = 0; y < l_s; y++) {
                    fread(&b, 1, 1, img);
                    fi[y + x * l_s].blue = b;
                    fread(&g, 1, 1, img);
                    fi[y + x * l_s].green = g;
                    fread(&r, 1, 1, img);
                    fi[y + x * l_s].red = r;
                }
                jump_s=3*(l_img-l_s)+padding_s;
                fseek(img,jump_s,SEEK_CUR);
            }
            corr=correlation(fi,s,h_s,l_s);
            if (corr > prag) {
                d[w].corelatie=corr;
                d[w].row=i;
                d[w].column=j;
                d[w].width=l_s;
                d[w].height=h_s;
                w++;

            }
        }
    }
    *n=w;
    return d;

}
int main()
{
    char nume_imagine_liniarizata[]="peppers.bmp";
    char nume_imagine_destinatie[]="peppers2.bmp";
    char nume_imagine_criptata[]="enc.bmp";
    char nume_imagine_decriptata[]="denc.bmp";
    char nume_cheie[]="secret_key.txt";
    struct pixel *a;
    unsigned int *h;
    incarca_liniarizata(&h, &a,nume_imagine_liniarizata);
    printf("Afisare test imagine initiala: \n");
    test_chi(nume_imagine_liniarizata,a);
    salveaza_liniarizata(h, a, nume_imagine_destinatie);
    criptare(nume_imagine_liniarizata, nume_imagine_criptata, nume_cheie,a);
    free(a);
    free(h);
    incarca_liniarizata(&h, &a,nume_imagine_criptata);
    printf("Afisare test imagine criptata: \n");
    test_chi(nume_imagine_criptata,a);
    decriptare(nume_imagine_decriptata,nume_imagine_criptata,nume_cheie);
    double prag=0.5;
    struct detectii *d;
    struct detectii *fi;
    unsigned int n,i,m=0;
    char test_initial[]="test.bmp";
    char test_border[]="test_new.bmp";
    char test_grey[]="test_grey.bmp";
    char cifra0[]="cifra0.bmp";
    char cifra1[]="cifra1.bmp";
    char cifra2[]="cifra2.bmp";
    char cifra3[]="cifra3.bmp";
    char cifra4[]="cifra4.bmp";
    char cifra5[]="cifra5.bmp";
    char cifra6[]="cifra6.bmp";
    char cifra7[]="cifra7.bmp";
    char cifra8[]="cifra8.bmp";
    char cifra9[]="cifra9.bmp";
    grayscale(test_initial,test_grey);
    incarca_liniarizata(&h,&a,test_initial);
    salveaza_liniarizata(h,a,test_border);
    fi=match(test_grey,cifra0,prag,&n);
    d=malloc(n*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra1,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra2,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra1,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra3,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra4,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra5,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra6,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra7,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra8,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    fi=match(test_grey,cifra9,prag,&n);
    d=realloc(d,(n+m)*sizeof(struct detectii));
    for(i=0;i<n;i++)
    {
        d[i+m]=fi[i];

    }
    free(fi);
    m=m+n;
    return 0;
}