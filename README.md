# FFTChordDetector
Project build around the MAX/MSP API for harmonic detection 
 with the help of Graham Wakefield wrapper for C++
Accelerate Framework for fourier analysis is used. 
I've made an FFT class for further reuse in platform that support the Accelerate Framework









Rapport du projet OOP/CPP
détecteur d’Harmonie


Contextualisation du projet.
Le but du projet est de montrer qu’il est possible de développer un outil analyseur d’harmonie
qui permet d’afficher l’accord en temps réels pour l’ environnement MAX for Live afin de permettre à l’artiste de savoir 
quelle gamme il compose sur des synthétiseur où le langage fréquentiel diffère de celui du langage musical.

Présentation du framework Accelerate
Le choix du framework Accelerate de Apple a été choisi pour l’algorithme de transformé rapide de fourier. 
La classe FFT pourra ainsi servir par la suite pour des applications iOS
#include <Accelerate/Accelerate.h>


Présentation du software MAX MSP
Max MSP est une platforme de développement de programmation graphique orienté sur la conception de synthèse sonore via MSP et des contrôles de ses paramètre via MAX. Tout les objets de la librairie sont développé en C. Le principe repose sur l’imbrication de plusieurs objets.  Il est donc possible de développer ses propres objets via des langages de programmation tiers tel que le Java-Script, le Java et le langage C. 

#include "ext.h"
#include "ext_obex.h"
#include "ext_common.h"
#include "commonsyms.h"
#include "z_dsp.h"
Présentation de la class FFT
    COMPLEX_SPLIT A;	// définit le complex buffer
    FFTSetup fftSetup;
    vDSP_Length log2n;
    uint32_t L;			//	taille du tableau 

Le Constructeur
    FFT(int L, float sr):L(L), sr(sr){
        
        log2n = log2(L);
        bw = (double) (2.0f / L) * (sr / 2.0f);
        
        
        fftSetup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
        
        A.realp = (float *) malloc(sizeof(float)* L/2);
        A.imagp = (float *) malloc(sizeof(float)* L/2);
        

       }   



La fonction forward 
(qui converti le signal temporel en domaine fréquentiel)
    
    void forward(float input[]){


        mag = new float[L/2];
        phase = new float[L/2];
        
        vDSP_vmul(input, 1, m_hammingWindow, 1, input, 1, L);

        vDSP_ctoz((COMPLEX*)input,2,&A,1,L/2);
        
        vDSP_fft_zrip(fftSetup, &A, 1, log2n, FFT_FORWARD);
    
        mag[0] = sqrt(A.realp[0]*A.realp[0]);

        //get phase
        vDSP_zvphas (&A, 1, phase, 1, L/2);
        phase[0] = 0;
        
        for(int i=1;i<L/2;i++)
            mag[i]=sqrt(A.realp[i]*A.realp[i]+A.imagp[i]*A.imagp[i]);
    
    
    }







La fonction frequencyToIndex
(qui permet de retrouver un indice à partir d’une fréquence)

    int frequencyToIndex(float freq){
        
        if(freq < bw/2)return 0;
        if(freq > sr/2 - bw/2)return L/2;
        
        float ratio =  freq / sr;
        int bin = round(L * ratio);
        return bin;
    }

La fonction ftom  

(qui convertit une fréquence en une Note MIDI.)

static inline int ftom (float f)
{
    return (int)(0.5 + 69.0 + 12.0 / logf (2.0) * logf (f / 440.0));
}

La fonction average_power

(qui correspond a un RMS de chaque bande)

float average_power(double * p, int start_bin, int stop_bin)
{
    float av = 0.0;
    for (int i = start_bin; i < stop_bin; i++)
    {
        av += p[i];
    }
    av /= (double)(stop_bin - start_bin);
    return av;
    }


Description du template maxcpp.h

Lorsque le constructeur de la class principale est appelé, il crée un object MAX/MSP qui réside a un pointeur vide via la fonction object_alloc




	static t_class * makeMaxClass(const char * classname) {
		common_symbols_init();
		t_class * c = class_new(classname, (method)MspCpp6<T>::maxcpp_create, (method)MspCpp6<T>::maxcpp_destroy, sizeof(T), NULL, A_GIMME, 0);
		class_dspinit(c);
		
		class_addmethod(c, (method)MspCpp6<T>::maxcpp_dsp64, "dsp64", A_CANT, 0);
		
		class_register(CLASS_BOX, c);
		MaxCppBase<T>::m_class = c;
		return c;
	}



static void maxcpp_destroy(t_object * x) {
		dsp_free((t_pxobject *)x);
		((T *)x)->~T();
	}




La fonction maxcpp_perform64 
en c++ est un template de la fonction MAX MSP

	static void maxcpp_perform64(t_object *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam) {
        ((T *)x)->perform(ins, numins, outs, numouts, sampleframes);
    }


// for DSP
#define REGISTER_PERFORM(CLASS, METHOD) object_method( \
	dsp64, \
	gensym("dsp_add64"), \
	(t_object *)this, \
	MaxMethodPerform64<&CLASS::METHOD>::call,\
	0, \
	NULL);







Description de la class principale

	MyFFT(t_symbol * sym, long ac, t_atom * av) {
        
        fft = new FFT(bufferSize, sys_getsr());

        dsp_setup((t_pxobject *)this, 2);                      
         // set up DSP for the instance and create signal inlets
        
        
        // prevent recycling of inputs for outputs
        m_ob.z_misc = Z_NO_INPLACE;
        this->out1 = outlet_new((t_object *)this, "list");       
        // outlet for displaying the bands
        this->out2 = outlet_new((t_object *)this, "list");       
        // outlet for displaying the histogram

        this->out3 = outlet_new((t_object *)this, "list");       
        // outlet for the sorted notes

        post("Chord Detector By Zirafkend");
        
	}
	
	~MyFFT() {
        delete(fft);
        delete buffer;
        post("object freed");
	}

Description du DSP


 	int cnt = 0;
    int size = 0;
    t_float *buffer;

	void perform(double **ins, long numins, double **outs, long numouts, long sampleframes) {
		
        t_double *inL = ins[0];         // First inlet`
        t_double *inR = ins[1];         // Second inlet`

    
        
// Buffering du signal de 64 trame en 4096

        if (size == 0) {
            buffer = new t_float[bufferSize];
        }
        
        for (int i = size; i < sampleframes + size; i++) {
            buffer[i] = (t_float) inL[i-size] + inR[i-size] / 2;        
			/// Packing Stereo signal from the two inlets into the buffer performing the average
        }
        
        size += sampleframes;
        if(size == bufferSize){
            size = 0;

            fft->forward(buffer);


// Définition préalable d’une tessiture qui correspond au piano.

            const int notelow = 28; const int notehigh = 103;           
			/// TODO add this constants to MAX/MSP arguments

            int i0 = fft->frequencyToIndex(mtof(notelow));
            int i1 = fft->frequencyToIndex(mtof(notehigh));


// Réduction du spectre FFT de 4096 bandes en 128 bandes midi.

            
            /// Averaging spectrum
            
            float * averageMidi = new float[128];
            int * n = new int[128];
            
            for (int i = 0; i < 128; i++) {
                n[i] = 0;
                averageMidi[i] = 0.0f;
            }
            
            for(int i = i0; i < i1; i++){
                
                float f = fft->indexToFrequency(i);
                int midi = ftom(f);
                if(midi >= 0 && midi < 128){
                    
                    float bin = fft->getBand(i);
                    averageMidi[midi] += sqrt(bin);
                    n[midi]++;
                }
                
            }
            
            // do the Root Mean Square
            //http://en.wikipedia.org/wiki/Root_mean_square
            for (int midi = 0; midi < 128; midi++) {
                if(n[midi] > 0){
                    averageMidi[midi] = averageMidi[midi] / n[midi]; // average
                    averageMidi[midi] = averageMidi[midi] * averageMidi[midi]; // square
                    
                }
            }



Triage des notes midi par le biais d’un tableau de fréquence qui correspond à une map en C++

			t_atom out3[12];

            
            std::map <int, float> mymap;
            
            for (int i = 0; i < 12; i++) {
                    mymap.insert(std::pair<int,float>(i, harmonicHistogram[i]));
            }
            
            // sort the map
            vector<pair<int, float> > mapcopy(mymap.begin(), mymap.end());
            sort(mapcopy.begin(), mapcopy.end(), less_second<int, float>());
            
            
            // Trim the map
            //mapcopy.erase ( mapcopy.begin(), mapcopy.end()-chordLength);
            
        // Affichage du résultat du tableau dans l’outlet de l’objet

            for (int i = 0; i < mapcopy.size(); i++) {
                atom_setlong(out3+i,mapcopy[i].first);

            }


        Code qui correspond à un compareTo en Java.          
  
     template <typename T1, typename T2>
     struct less_second {
     typedef pair<T1, T2> type;
     bool operator ()(type const& a, type const& b) const {
     return a.second < b.second;
     }
     };
     





Référence


https://developer.apple.com/library/ios/documentation/Accelerate/Reference/vDSPRef/Reference/reference.html

http://stackoverflow.com/questions/18931905/using-apples-accelerate-framework-fft-hann-windowing-and-overlapping

https://code.google.com/p/maxcpp/

http://cycling74.com/products/max/

https://www.ableton.com/en/live/max-for-live/



