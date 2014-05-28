#include "maxcpp6.h"
#include "FFT.h"

// inherit from the MSP base class, template-specialized for myself:

#define BUFFER_SIZE 4096

static inline double mtof (int midi)
{
    return exp (log (440.0) + (double)(midi - 69) * log (2.0) / 12.0);
}

static inline int ftom (double f)
{
    return (int)(0.5 + 69.0 + 12.0 / log (2.0) * log (f / 440.0));
}


class MyFFT : public MspCpp6<MyFFT> {
public:
    
    FFT *fft;
    void * out1;                                    // Void pointer to an outlet


	MyFFT(t_symbol * sym, long ac, t_atom * av) {
		fft = new FFT(BUFFER_SIZE, sys_getsr());
        dsp_setup((t_pxobject *)this, 2);                       // set up DSP for the instance and create signal inlets
        // prevent recycling of inputs for outputs
        m_ob.z_misc = Z_NO_INPLACE;
        this->out1 = outlet_new((t_object *)this, "list");       // outlet for displaying the bands
        post("Chord Detector By Zirafkend");
	}
	
	~MyFFT() {
        delete(fft);
        delete buffer;
		post("object freed"); 
	}	
	
	// methods:
	void bang(long inlet) { 
		post("bang in inlet %i!", inlet); 
	}
	void test(long inlet, t_symbol * s, long ac, t_atom * av) { 
		post("%s in inlet %i (%i args)", s->s_name, inlet, ac);
	}
    
    /// Global variable for the buffering into the fft size
    
    int cnt = 0;
    int size = 0;
    t_float *buffer;
    
	
	void perform(double **ins, long numins, double **outs, long numouts, long sampleframes) {
		
        t_double *inL = ins[0];         // First inlet`
        t_double *inR = ins[1];         // Second inlet`
        
        
        
        if (size == 0) {
            buffer = new t_float[BUFFER_SIZE];
        }
        
        for (int i = size; i < sampleframes + size; i++) {
            buffer[i] = (t_float) inL[i-size] + inR[i-size] / 2;
            /// Packing Stereo signal from the two inlets into the buffer performing the average
        }
        
        size += sampleframes;
        if(size == BUFFER_SIZE){
            size = 0;
            
            fft->forward(buffer);
            
            
            // Cepstrum
            
            float *tempFFT = new float[BUFFER_SIZE];
            
            for (int i = 0; i < BUFFER_SIZE; i++) {
                tempFFT[i] = log(fft->get_band(i));
            }
            
            fft->forward(fft->inverse(tempFFT));
            
            
            /// Averaging spectrum
            
            
            //const int notelow = 28; const int notehigh = 103;           /// TODO add this constants to max arguments
            float * averageMidi = new float[128];
            int * n = new int[128];
            
            for (int i = 0; i < 128; i++) {
                n[i] = 0;
                averageMidi[i] = 0.0f;
            }
            
            
            int i0 = fft->frequencyToIndex(mtof(0));
            int i1 = fft->frequencyToIndex(mtof(127));
            
            
            
            for(int i = i0; i < i1; i++){
                
                float f = fft->indexToFrequency(i);
                int midi = ftom(f);
                if(midi >= 0 && midi < 128){
                    averageMidi[midi] += sqrt(fft->get_band(i));
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
            
            
            
            t_atom out1[128];
            
            for (int i=0; i < 128; i++) {
                atom_setlong(out1+i,averageMidi[i]);
                
            }
            outlet_list(this->out1,0L,128,out1);
            
            
            
            /*
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
             
             for (int i = 0; i < mapcopy.size(); i++) {
             atom_setlong(out3+i,mapcopy[i].first);
             
             }
             
             
             outlet_list(this->out3,0L,12,out3);
             
             
             
             
             template <typename T1, typename T2>
             struct less_second {
             typedef pair<T1, T2> type;
             bool operator ()(type const& a, type const& b) const {
             return a.second < b.second;
             }
             };
             
             */
            
            
        }
        
        
        
        
        
	}
	
//	// optional method: gets called when the dsp chain is modified
//	// if not given, the MspCpp will use Example::perform by default
//	void dsp(t_object * dsp64, short *count, double samplerate, long maxvectorsize, long flags) { 
//		// specify a perform method manually:
//		REGISTER_PERFORM(Example, perform);
//	}
};

C74_EXPORT int main(void) {
	// create a class with the given name:
	MyFFT::makeMaxClass("example~");
	REGISTER_METHOD(MyFFT, bang);
	REGISTER_METHOD_GIMME(MyFFT, test);
}