//
//  FFT.h
//  simplemsp~
//
//  Created by Nikola Markovic on 20/12/13.
//
//

#ifndef __simplemsp___FFT__
#define __simplemsp___FFT__

#include <iostream>
#include <Accelerate/Accelerate.h>




int isPowerOfTwo (unsigned long x) { return ((x != 0) && !(x & (x - 1))); }



class FFT{
    
protected:
    
    float sr;
    uint32_t N, N_2;
    
    // Define complex buffer
    COMPLEX_SPLIT A;
    COMPLEX_SPLIT B;
    
    FFTSetup fftSetup;

    float * m_hammingWindow;
    float *mag;
    float *phase;
    
    float bw;
    
    vDSP_Length log2n;
    
public:
    

    
    FFT(int N, float sr):N(N), sr(sr){
        
        log2n = log2(N);
        
        N_2 = N / 2;
        
        bw =  (2.0f / N) * (sr / 2.0f);
        
        
        fftSetup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
        
        m_hammingWindow = (float *) malloc(sizeof(float) * N);
        
        vDSP_hamm_window(m_hammingWindow, N, vDSP_HANN_NORM);
        
        
        A.realp = (float *) malloc(sizeof(float)* N_2);
        A.imagp = (float *) malloc(sizeof(float)* N_2);
        
        B.realp = (float*) malloc(sizeof(float) * N_2);
        B.imagp = (float*) malloc(sizeof(float) * N_2);
        
    }
    
    ~FFT(){
        delete A.realp;
        delete A.imagp;
        
        delete B.realp;
        delete B.imagp;
        vDSP_destroy_fftsetup(fftSetup);

    }
    
    uint32_t specSize(){
        return N;
    }
    
    float get_band(int i){
        if(i < 0) i = 0;
        if(i > N - 1 ) i = N -1;
        return mag[i];
    }
    
   
    
    void forward(float input[]){


        mag = new float[N_2];
        phase = new float[N_2];
        
        vDSP_vmul(input, 1, m_hammingWindow, 1, input, 1, N);

        vDSP_ctoz((COMPLEX*)input,2,&A,1,N_2);
        
        vDSP_fft_zrip(fftSetup, &A, 1, log2n, FFT_FORWARD);
    
        mag[0] = sqrt(A.realp[0]*A.realp[0]);

        //get phase
        vDSP_zvphas (&A, 1, phase, 1, N_2);
        phase[0] = 0;
        
        for(int i=1;i<N_2;i++)
            mag[i]=sqrt(A.realp[i]*A.realp[i]+A.imagp[i]*A.imagp[i]);
    
    
    }
    
    
    
    float* inverse(float buffer[]){
        
        float *output = new float[N];

        //after done with possible phase and mag processing re-pack the vectors in VDSP format
        B.realp[0] = mag[0];
        B.imagp[0] = mag[N/2 - 1];;
        
        //unwrap, process & re-wrap phase
        for(int i = 1; i < N_2; i++){

            phase[i] -= 2*M_PI*i * sr/N;
            phase[i] -= M_PI / 2 ;
            phase[i] += 2*M_PI*i * sr/N;
        }
        
        //construct real & imaginary part of the output packed vector (input to ifft)
        for(int i = 1; i < N_2; i++){
            B.realp[i] = mag[i] * cosf(phase[i]);
            B.imagp[i] = mag[i] * sinf(phase[i]);
        }
        
        vDSP_fft_zrip(fftSetup, &B, 1, log2n, FFT_INVERSE);
        //scale factor
        float scale = 1.0 / (2*N);
        
        //scale values
        vDSP_vsmul(B.realp, 1, &scale, B.realp, 1, N_2);
        vDSP_vsmul(B.imagp, 1, &scale, B.imagp, 1, N_2);
        
        //unpack B to real interleaved output
        vDSP_ztoc(&B, 1, (COMPLEX *) output, 2, N_2);
        
        return output;

    }
    

    int frequencyToIndex(float freq){
        
        if(freq < bw/2)return 0;
        if(freq > sr/2 - bw/2)return N_2;
        
        double ratio =  freq / sr;
        int bin = round(N * ratio);
        return bin;
    }
    
    
    float indexToFrequency(int i){
        if ( i == 0 ) return bw * 0.25f;
        if ( i == (N_2)){
            float lastBinBeginFreq = (sr / 2) - (bw / 2);
            float binHalfWidth = bw * 0.25f;
            return lastBinBeginFreq + binHalfWidth;
        }
        return i*bw;
        
    }

    int maxBin(int start_bin, int stop_bin){
        int max = -1;
        int imax;

        for (int i = start_bin; i < stop_bin; i++)
        {
            if (mag[i] > max)
            {
                max = mag[i];
                imax = i;
            }
        }
        if(imax == -1) return -1;
        return max;
    }

    float average_power(int start_bin, int stop_bin)
    {
        float av = 0.0;
        for (int i = start_bin; i < stop_bin; i++)
        {
            av += mag[i];
        }
        av /= (double)(stop_bin - start_bin);
        return av;

    }



};


void phase_correction(int buffer_size, double *p,  double *p0 ,double *dphi ,double *ph0,double *ph1, int hop)
{
    for (int i = 0; i < (buffer_size/2+1); ++i) // full span
    {
        double twopi = 2.0 * M_PI;
        //double dphi;
        dphi[i] = ph1[i] - ph0[i]
                - twopi * (double)i / (double)buffer_size * (double)hop;
        for (; dphi[i] >= M_PI; dphi[i] -= twopi);
        for (; dphi[i] < -M_PI; dphi[i] += twopi);

        // frequency correction
        // NOTE: freq is (i / len + dphi) * samplerate [Hz]
        dphi[i] = dphi[i] / twopi / (double)hop;

        // backup the phase for the next step
        p0  [i] = p   [i];
        ph0 [i] = ph1 [i];

        // then, average the power for the analysis
        p[i] = 0.5 *(sqrt (p[i]) + sqrt (p0[i]));
        p[i] = p[i] * p[i];
    }
}




class BeatDetector : FFT {
    
    
    float *magnitude_average;
    float *averageEnergy;
    float *fftSubbands;
    float *beatValueArray;
    float *fftVariance;
    float** energyHistory;
    
    int SUBBANDS = 128;
    int ENERGY_HISTORY = 64;
    
    
    BeatDetector(int N, float sr):FFT(N, sr){
        magnitude_average = new float[N];
        
        for(int i = 0; i < SUBBANDS; i++)
        {
            for(int l = 0; l < ENERGY_HISTORY; l++){
                energyHistory[i][l] = 0;
            }
            fftSubbands[i] = 0;
            averageEnergy[i] = 0;
            fftVariance[i] = 0;
            beatValueArray[i] = 0;
        }
    
    };
    

    
    
    void computeRMS(float *input){
        forward(input);
        for (int i = 0; i < N; i++) {
            float x = 0.085;
            magnitude_average[i] = (mag[i] * x) + (mag[i] * (1 - x));
        }
    }
    
    
    void computeSubands(){
        
        
        for(int i = 0; i < SUBBANDS; i++) {
            fftSubbands[i] = 0;
            
            for(int b = 0; b < N/SUBBANDS; b++) {
                fftSubbands[i] +=  get_band(i*(N/SUBBANDS)+b);
            }
            fftSubbands[i] = fftSubbands[i]*(float)SUBBANDS/(float)N;
            
            for(int b=0; b < N/SUBBANDS; b++) {
                fftVariance[i] +=  powf(get_band(i*(N/SUBBANDS)+b)- fftSubbands[i], 2);
            }
            fftVariance[i] = fftVariance[i]*(float)SUBBANDS/(float)N;
            
            beatValueArray[i] = (-0.0025714*fftVariance[i])+1.35;
        }
        
        for(int i = 0; i < SUBBANDS; i++) {
            averageEnergy[i] = 0;
            for(int h = 0; h < ENERGY_HISTORY; h++) {
                averageEnergy[i] += energyHistory[i][h];
            }
            averageEnergy[i] /= ENERGY_HISTORY;
        }

    
    }

    bool isBeat(int subband)
    {
        return fftSubbands[subband] > averageEnergy[subband]*beatValueArray[subband];
    }
    
    
    bool isBeatRange(int low, int high, int threshold)
    {
        int num = 0;
        for(int i = low; i < high+1; i++) {
            if(isBeat(i)) {
                num++;
            }
        }
        return num > threshold;
    }
    
    
    
    bool isKick()
    {
        return isBeat(0);
    }
    
    bool isSnare()
    {
        int low = 1;
        int hi = SUBBANDS/3;
        int thresh = (hi-low)/3;
        return isBeatRange(low, hi, thresh);
    }
    
    bool isHat()
    {
        int low = SUBBANDS/2;
        int hi = SUBBANDS-1;
        int thresh = (hi-low)/3;
        return isBeatRange(low, hi, thresh);
    }
    
    
    
};



/*



*/

#endif /* defined(__simplemsp___FFT__) */
