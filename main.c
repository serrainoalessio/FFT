#include <stdio.h> // for printf
#include <stdlib.h> // for malloc
#include <string.h> // for memcpy
#include <math.h> // for sin & cos

#define PI 3.14159265358979323 // obviously pi has infinite decimals
#define PI_DOUBLE PI*2

typedef struct {
    double real, imag;
} complex_t;

double cmplx_real(complex_t * op) {
    return op->real;
}

double cmplx_imag(complex_t * op) {
    return op->imag;
}

void cmplx_add(complex_t * res, complex_t * op1, complex_t * op2) { // sum
    res->real = op1->real + op2->real;
    res->imag = op1->imag + op2->imag;
}

void cmplx_sub(complex_t * res, complex_t * op1, complex_t * op2) { // difference
    res->real = op1->real - op2->real;
    res->imag = op1->imag - op2->imag;
}

void cmplx_mul(complex_t * res, complex_t * op1, complex_t * op2) { // product
    double tmp; // this to handle the case res == one of the ops
    tmp       = (op1->real)*(op2->imag) + (op1->imag)*(op2->real);
    res->real = (op1->real)*(op2->real) - (op1->imag)*(op2->imag);
    res->imag = tmp;
}

void cmplx_div(complex_t *res, complex_t * op1, complex_t * op2) {
    double den = pow(op2->real, 2) + pow(op2->imag, 2);
    double temp = (op1->imag)*(op2->real) - (op1->real)*(op2->imag);
    res->real   = (op1->real)*(op2->real) + (op1->imag)*(op2->imag);
    res->imag = temp;
    res->real /= den;
    res->imag /= den;
}

void cmplx_cpy(complex_t * dest, complex_t * source) { // copy
    memcpy(dest, source, sizeof*dest);
}

void cmplx_cnj(complex_t * res, complex_t * op) { // complex conjugate
    res->real = +(op->real);   // the complex conjugate of  A + iB is
    res->imag = -(op->imag);   //                          A - iB
}

void cmplx_swp(complex_t * res, complex_t * op) { // swaps real and imag
    double swap = (op->real);
    res->real   = (op->imag);
    res->imag   = swap;
}

int invert(int n, int len) {
    int i, res = 0;
    for (i = 0; i < len; i++)
        if (n & (1 << i))
            res |= (1 << len - (i + 1));
    return res;
}

void reorder(complex_t * datablock, int len) { // reorders data for fft
    int i, j, bits;
    complex_t swap;
    bits = __builtin_ctz(len); // counts ending zeros
                               // len must be a power of two, so it is log2(len)
    for (i = 0; i < len; i++) {
        j = invert(i, bits);
        if (i < j) {
            swap = datablock[i];
            datablock[i] = datablock[j];
            datablock[j] = swap;
        }
    }
}

void fft(complex_t * datablock, int len) { // recursively computes fft on an already ordered block
    int k, m = len/4;
    complex_t twiddle, temp;
    double theta;

    if (len <= 1) { // FFT of one lenght is defined as the input, so no changes!
        return;
    } else if (len == 2) {
        cmplx_cpy(&twiddle, datablock + 1); // saves temporany data
        cmplx_cpy(&temp   , datablock + 0); // saves temporany data
        cmplx_add(datablock + 0, &temp, &twiddle);
        cmplx_sub(datablock + 1, &temp, &twiddle);
        return;
    } /* else if (len == 4) {
        cmplx_cpy(&twiddle, datablock + 1); // saves temporany data
        cmplx_cpy(&temp   , datablock + 0); // saves temporany data
        cmplx_add(datablock + 0, &temp, &twiddle);
        cmplx_sub(datablock + 1, &temp, &twiddle);
        return;
    } */

    // split four radix fft
    fft(datablock + 0*m, m); // fft of 0s (mod 4)
    fft(datablock + 1*m, m); // fft of 1s (mod 4)
    fft(datablock + 2*m, m); // fft of 2s (mod 4)
    fft(datablock + 3*m, m); // fft of 3s (mod 4)

    for (k = 0; k < m; k++) {
        theta = -PI_DOUBLE*(k+0*m)/len; // now computes the angle and the twiddle factor
        twiddle.real = cos(theta);
        twiddle.imag = sin(theta);
        cmplx_mul(&twiddle, &twiddle, datablock + k + 1*m); // twiddle computed!
        cmplx_cpy(&temp, datablock + k + 0*m); // saves temporany data
        cmplx_add(datablock + k + 0*m, &temp, &twiddle);

        theta = -PI_DOUBLE*(k+1*m)/len; // now computes the angle and the twiddle factor
        twiddle.real = cos(theta);
        twiddle.imag = sin(theta);
        cmplx_mul(&twiddle, &twiddle, datablock + k + 1*m); // twiddle computed!
        cmplx_cpy(&temp, datablock + k + 2*m); // saves temporany data
        cmplx_add(datablock + k + 1*m, &temp, &twiddle);

        theta = -PI_DOUBLE*(k+2*m)/len; // now computes the angle and the twiddle factor
        twiddle.real = cos(theta);
        twiddle.imag = sin(theta);
        cmplx_mul(&twiddle, &twiddle, datablock + k + 3*m); // twiddle computed!
        cmplx_cpy(&temp, datablock + k + 0*m); // saves temporany data
        cmplx_add(datablock + k + 2*m, &temp, &twiddle);

        theta = -PI_DOUBLE*(k+3*m)/len; // now computes the angle and the twiddle factor
        twiddle.real = cos(theta);
        twiddle.imag = sin(theta);
        cmplx_mul(&twiddle, &twiddle, datablock + k + 3*m); // twiddle computed!
        cmplx_cpy(&temp, datablock + k + 2*m); // saves temporany data
        cmplx_add(datablock + k + 3*m, &temp, &twiddle);
    }
}

void ifft(complex_t * datablock, int len) { // inverse fft
    int i;
    complex_t scale = (complex_t){.real = len, .imag = 0};

    //first calculates the complex conjugate of every input
    for (i = 0; i < len; i++)
        cmplx_cnj(datablock + i, datablock + i);
    //computes a normal forward fft
    fft(datablock, len);

    //computes the complex conjugate of the result
    for (i = 0; i < len; i++)
        cmplx_cnj(datablock + i, datablock + i);

    //scale down the result
    for (i = 0; i < len; i++)
        cmplx_div(datablock + i, datablock + i, &scale);
}

void window(complex_t * datablock, double * win, int len) { // windowing function
    int i;
    complex_t temp;

    for (i = 0; i < len; i++) {
        temp = (complex_t){.real = win[i], .imag = 0.0};
        cmplx_mul(datablock + i, datablock + i, &temp);
    }
}

void unwindow(complex_t * datablock, double * win, int len) { // unuseful
    int i;
    complex_t temp;

    for (i = 0; i < len; i++) {
        temp = (complex_t){.real = win[i], .imag = 0.0};
        cmplx_div(datablock + i, datablock + i, &temp);
    }
}

#define LEN 150
#define BAND 16
#define ZERO 8

#define DISPLAY_MODULO 0
#define DISPLAY_REAL   1
#define DISPLAY_IMAG   2

#define DISPLAY_HALF_MODULO 3

void draw(double * datablock, int len) { // draws the funcion on the screen
    int band, i, j, max, min;
    double idx0, idx1, threeshold_min,
                       threeshold_max,
                       lband, data;

    for (i = 1, max = min = 0; i < len; i++) {
        if (datablock[max] < datablock[i])
            max = i;
        if (datablock[min] > datablock[i])
            min = i;
    }

    lband = (datablock[max] - datablock[min])/BAND;

    for (band = BAND+1; band >= 0; band--) {
        threeshold_min = datablock[min]+band*lband;
        threeshold_max = datablock[min]+(band+1)*lband;
        for (i = 0; i < LEN; i++) {
            idx0 = ceil(len*i/LEN);
            idx1 = floor(len*(i+1)/LEN);
            data = 0.0;
            for (j = idx0; j < idx1; j++)
                data += datablock[j];
            data /= (idx1 - idx0);

            if ((data > threeshold_min) && (data <= threeshold_max)) {
                printf("*");
            } else {
                if (band == ZERO)
                    printf("-");
                else
                    printf(" ");
            }
        }
        printf("\n");
    }
}

void cmplx_draw(complex_t * datablock, int len, int mode) { // draws the modulo of the datablock
    double * data;
    int i;

    if (mode == DISPLAY_HALF_MODULO)
        data = malloc((len/2)*sizeof*data);
    else
        data = malloc(len*sizeof*data);

    for (i = 0; i < len; i++) {
        switch (mode) {
          case DISPLAY_MODULO:
            data[i] = sqrt(pow(cmplx_real(datablock + i), 2) + pow(cmplx_imag(datablock + i), 2));
            break;
          case DISPLAY_REAL:
            data[i] = cmplx_real(datablock + i);
            break;
          case DISPLAY_IMAG:
            data[i] = cmplx_imag(datablock + i);
            break;
        case DISPLAY_HALF_MODULO:
            if (i < len/2)
                data[i] = sqrt(pow(cmplx_real(datablock + i), 2) + pow(cmplx_imag(datablock + i), 2));
            break;
          default:
            data[i] = 0;
            break; // this should be an error
        }
    }

    if (mode == DISPLAY_HALF_MODULO)
        draw(data, len/2);
    else
        draw(data, len);

    free(data);
}

int main(int argc, char * argv) {
    complex_t * data; // data
    double * win; // windowing function
    int i, len;

    do {
        printf("insert the lenght: ");
        scanf("%d", &len);
    } while (len & (len - 1)); // checks that len is a power of two

    data = malloc(len*sizeof*data);
    win  = malloc(len*sizeof*win);

    // now inits datas:
    for (i = 0; i < len; i++) {
        /* printf("Insert real & imag [%d]: ", i);
         * scanf("%lf", &data[i].real);
         * scanf("%lf", &data[i].imag);
         */

        data[i].real = sin(10.0*PI_DOUBLE*i/len) + sin(28.0*PI_DOUBLE*i/len);
        data[i].imag = 0;
    }

    // now init window
    for (i = 0; i < len; i++) { // Hamming
        win[i] = 0.53836 - 0.46164*cos(PI_DOUBLE*i/(len-1));
        win[i] = 1; // no windowing
    }

    printf("\n");


    printf("input:\n");
    cmplx_draw(data, len, DISPLAY_REAL);
/*
    for (i = 0; i < len; i++)
        printf("%d: (%lf; %lf) %lf\n", i, data[i].real, data[i].imag, pow(data[i].real, 2) + pow(data[i].imag, 2) ); */
    printf("\n");

    window(data, win, len); // windows the function
    reorder(data, len); // reorders data
    fft(data, len); // computes forward fft

    printf("output:\n");
    cmplx_draw(data, len, DISPLAY_HALF_MODULO);
/*
    for (i = 0; i < len; i++)
        printf("%d: (%lf; %lf) %lf\n", i, data[i].real, data[i].imag, pow(data[i].real, 2) + pow(data[i].imag, 2) ); */
    printf("\n");

    reorder(data, len); // twice reorder execution generates the starting data
    ifft(data, len); // inverse fft
    unwindow(data, win, len); // <--- this may lose mutch precision

    printf("check:\n");
    cmplx_draw(data, len, DISPLAY_REAL);
/*
    printf("check:\n");
    for (i = 0; i < len; i++)
        printf("%d: (%lf; %lf) %lf\n", i, data[i].real, data[i].imag, pow(data[i].real, 2) + pow(data[i].imag, 2) );
*/

    free(data);
    free(win);

    return 0; // no errors
}
