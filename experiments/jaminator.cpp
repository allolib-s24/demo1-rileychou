// Copied from 08_SubSyn.cpp

#include <cstdio>  // for printing to stdout

#include <vector>
#include <random>
#include <iostream>
#include <map>

#include "Gamma/Analysis.h"
#include "Gamma/Effects.h"
#include "Gamma/Envelope.h"
#include "Gamma/Gamma.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Types.h"
#include "Gamma/Spatial.h"
#include "Gamma/SamplePlayer.h"
#include "Gamma/DFT.h"

#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/scene/al_PolySynth.hpp"
#include "al/scene/al_SynthSequencer.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"

#include "al/app/al_App.hpp"
#include "al/io/al_Imgui.hpp"

using namespace gam;
using namespace al;
using namespace std;
#define FFT_SIZE 4048

class Snare : public SynthVoice {
 public:
  // Unit generators
  gam::NoiseWhite<> mNoise; // Use white noise for drum sound
  gam::Env<1> mAmpEnv;       // Use a simple envelope with only one segment
  
  // Initialize voice
  void init() override {
    // Intialize envelope
    mAmpEnv.curve(0); // linear curve
    mAmpEnv.levels(1, 0); // start with full amplitude and then decay to 0 immediately
    
    // Create parameters
    createInternalTriggerParameter("amplitude", 0.8, 0.0, 1.0);
    createInternalTriggerParameter("decayTime", 0.1, 0.01, 2.0);
  }
  //look at subtractive synthesis examples -- can start with this,
  // filters - filter out high and low frequencies
  // Audio processing function
  void onProcess(AudioIOData& io) override {
    while (io()) {
      float s1 = mNoise() * mAmpEnv() * getInternalParameterValue("amplitude");
      io.out(0) += s1;
      io.out(1) += s1;
    }
    
    // Free the voice when the envelope is done
    if (mAmpEnv.done()) free();
  }

  // Triggering functions
  void onTriggerOn() override {
    // Set the decay time when the note is triggered
    mAmpEnv.lengths()[0] = getInternalParameterValue("decayTime");
    mAmpEnv.reset();
  }

  void onTriggerOff() override {}
};


class BassDrum : public SynthVoice {
 public:
  // Unit generators
  gam::Sine<> mOsc;           // Sine wave oscillator for the bass drum tone
  gam::Env<1> mAmpEnv;         // Envelope for controlling the amplitude

  // Initialize voice
  void init() override {
    // Initialize envelope
    mAmpEnv.curve(0);          // linear curve
    mAmpEnv.levels(1, 0);      // start with full amplitude and then decay to 0 immediately
    
    // Create parameters
    createInternalTriggerParameter("frequency", 80.0, 20.0, 25.0); // Set the initial frequency to 80 Hz
    createInternalTriggerParameter("amplitude", 0.8, 0.0, 1.0);
    createInternalTriggerParameter("decayTime", 0.1, 0.01, 2.0);
  }

  // Audio processing function
  void onProcess(AudioIOData& io) override {
    float f = getInternalParameterValue("frequency");
    mOsc.freq(f);
    while (io()) {
      float s1 = mOsc() * mAmpEnv() * getInternalParameterValue("amplitude");
      io.out(0) += s1;
      io.out(1) += s1;
    }
    
    // Free the voice when the envelope is done
    if (mAmpEnv.done()) free();
  }

  // Triggering functions
  void onTriggerOn() override {
    // Set the frequency and decay time when the note is triggered
    // mOsc.freq(getInternalParameterValue("frequency"));
    mAmpEnv.lengths()[0] = getInternalParameterValue("decayTime");
    mAmpEnv.reset();
  }

  void onTriggerOff() override {}
};

class Cymbal : public SynthVoice {
 public:
  // Unit generators
  gam::NoiseWhite<> mNoise; // White noise generator for cymbal sound
  gam::Biquad<> mFilter;    // Biquad filter for shaping the cymbal sound
  gam::Env<3> mAmpEnv;       // Envelope for controlling the amplitude

  // Initialize voice
  void init() override {
    // Initialize filter
    mFilter.type(gam::HIGH_PASS);
    mFilter.freq(1000); // Set initial cutoff frequency
    mFilter.res(2.0);   // Set initial resonance

    // Initialize envelope
    mAmpEnv.curve(0);  // linear curve
    mAmpEnv.levels(0, 1, 0.2); // start with no amplitude, peak quickly, then decay
    
    // Create parameters
    createInternalTriggerParameter("cutoffFreq", 1000.0, 20.0, 5000.0); // Set the initial cutoff frequency
    createInternalTriggerParameter("resonance", 2.0, 0.1, 10.0);        // Set the initial resonance
    createInternalTriggerParameter("amplitude", 0.8, 0.0, 1.0);
    createInternalTriggerParameter("attackTime", 0.001, 0.01, 1.0);
    createInternalTriggerParameter("decayTime", 0.1, 0.01, 5.0);
  }

  // Audio processing function
  void onProcess(AudioIOData& io) override {
    while (io()) {
      // Generate noise
      float s1 = mNoise() * mAmpEnv() * getInternalParameterValue("amplitude");
      
      // Filter the noise
      mFilter.freq(getInternalParameterValue("cutoffFreq"));
      mFilter.res(getInternalParameterValue("resonance"));
      s1 = mFilter(s1);
      
      // Write to output channels
      io.out(0) += s1;
      io.out(1) += s1;
    }
    
    // Free the voice when the envelope is done
    if (mAmpEnv.done()) free();
  }

  // Triggering functions
  void onTriggerOn() override {
    // Set envelope parameters when note is triggered
    mAmpEnv.lengths()[0] = getInternalParameterValue("attackTime");
    mAmpEnv.lengths()[2] = getInternalParameterValue("decayTime");
    mAmpEnv.reset();
  }

  void onTriggerOff() override {}
};


class DemoVoice : public SynthVoice {
  public:
  gam::Pan<> mPan;
  gam::Sine<> mOsc;
  gam::Env<3> mAmpEnv;
  gam::EnvFollow<> mEnvFollow;

  gam::STFT stft = gam::STFT(FFT_SIZE, FFT_SIZE/4, 0, gam::HANN, gam::MAG_FREQ);

  Mesh mDemoVoice;
  vector<float> spectrum;
  Mesh mMesh;

  void init() override {
    mAmpEnv.curve(0); // linear segments
    mAmpEnv.levels(0,1,1,0);
    mAmpEnv.sustainPoint(2); // Make point 2 sustain until a release is issued

    // Declare the size of the spectrum
    spectrum.resize(FFT_SIZE);
    mDemoVoice.primitive(Mesh::LINE_LOOP);
    // mDemoVoice.primitive(Mesh::POINTS);

    addDisc(mMesh, 1.0, 30);

    createInternalTriggerParameter("amplitude", 0.3, 0.0, 1.0);
    createInternalTriggerParameter("frequency", 60, 20, 5000);
    createInternalTriggerParameter("attackTime", 0.1, 0.01, 3.0);
    createInternalTriggerParameter("releaseTime", 0.1, 0.1, 10.0);
    createInternalTriggerParameter("pan", 0.0, -1.0, 1.0);
  }

  void onProcess(AudioIOData& io) override {
    // Get the values from the parameters and apply them to the corresponding
    // unit generators. You could place these lines in the onTrigger() function,
    // but placing them here allows for realtime prototyping on a running
    // voice, rather than having to trigger a new voice to hear the changes.
    // Parameters will update values once per audio callback because they
    // are outside the sample processing loop.
    float f = getInternalParameterValue("frequency");
    mOsc.freq(f);
    mAmpEnv.lengths()[0] = getInternalParameterValue("attackTime");
    mAmpEnv.lengths()[2] = getInternalParameterValue("releaseTime");
    mPan.pos(getInternalParameterValue("pan"));
    while(io()){
      float s1 = mOsc() * mAmpEnv() * getInternalParameterValue("amplitude");
      float s2;
      mPan(s1, s1, s2);
      mEnvFollow(s1);
      io.out(0) += s1;
      io.out(1) += s2;

      if(stft(s1)){
        // loop through all frequency bins and scale the complex sample
        for(unsigned k = 0; k < stft.numBins(); ++k){
          spectrum[k] = tanh(pow(stft.bin(k).real(), 1.3));
        }
      }
    }
    // We need to let the synth know that this voice is done
    // by calling the free(). This takes the voice out of the
    // rendering chain
    if(mAmpEnv.done()) free();
  }

  void onProcess(Graphics& g) override {
    // Figure out graphics for chords later

    float frequency = getInternalParameterValue("frequency");
    float amplitude = getInternalParameterValue("amplitude");

    mDemoVoice.reset();

    for(int i = 0; i < FFT_SIZE / 90; i++){
      mDemoVoice.color(HSV(frequency/500 - spectrum[i] * 50));
      mDemoVoice.vertex(i, spectrum[i], 0.0);
    }

    for(int i = -5; i <= 4; i++) {
      g.meshColor();
      g.pushMatrix();
      // g.translate(-0.5f, 1, -10);
      // g.translate(cos(static_cast<double>(frequency)), sin(static_cast<double>(frequency)), -4);
      g.translate(i, -2.7, -10);
      g.scale(50.0/FFT_SIZE, 250, 1.0);
      // g.pointSize(1 + 5 * mEnvFollow.value() * 10);
      g.lineWidth(1 + 5 * mEnvFollow.value() * 100);
      g.draw(mDemoVoice);
      g.popMatrix();
    }
  }

  void onTriggerOn() override {
    mAmpEnv.reset();
  }

  void onTriggerOff() override {
    mAmpEnv.release();
  }
};

class PluckedString : public SynthVoice
{
public:
    float mAmp;
    float mDur;
    float mPanRise;
    gam::Pan<> mPan;
    gam::NoiseWhite<> noise;
    gam::Decay<> env;
    gam::MovingAvg<> fil{2};
    gam::Delay<float, gam::ipl::Trunc> delay;
    gam::ADSR<> mAmpEnv;
    gam::EnvFollow<> mEnvFollow;
    gam::Env<2> mPanEnv;
    gam::STFT stft = gam::STFT(FFT_SIZE, FFT_SIZE / 4, 0, gam::HANN, gam::MAG_FREQ);
    gam::Chorus<> chorus;
    // This time, let's use spectrograms for each notes as the visual components.
    Mesh mSpectrogram;
    vector<float> spectrum;
    double a = 0;
    double b = 0;
    double timepose = 10;
    // Additional members
    Mesh mMesh;

    virtual void init() override
    {
        // Declare the size of the spectrum
        spectrum.resize(FFT_SIZE / 2 + 1);
        // mSpectrogram.primitive(Mesh::POINTS);
        mSpectrogram.primitive(Mesh::LINE_STRIP);
        mAmpEnv.levels(0, 1, 1, 0);
        mPanEnv.curve(4);
        env.decay(0.1);
        delay.maxDelay(1. / 27.5);
        delay.delay(1. / 440.0);

        addDisc(mMesh, 1.0, 30);
        createInternalTriggerParameter("amplitude", 0.1, 0.0, 1.0);
        createInternalTriggerParameter("frequency", 60, 20, 5000);
        createInternalTriggerParameter("attackTime", 0.001, 0.001, 1.0);
        createInternalTriggerParameter("releaseTime", 3.0, 0.1, 10.0);
        createInternalTriggerParameter("sustain", 0.7, 0.0, 1.0);
        createInternalTriggerParameter("Pan1", 0.0, -1.0, 1.0);
        createInternalTriggerParameter("Pan2", 0.0, -1.0, 1.0);
        createInternalTriggerParameter("PanRise", 0.0, 0, 3.0); // range check
    }

    //    void reset(){ env.reset(); }

    float operator()()
    {
        return (*this)(noise() * env());
    }
    float operator()(float in)
    {
        return delay(
            fil(delay() + in));
    }

    virtual void onProcess(AudioIOData &io) override
    {

        while (io())
        {
            mPan.pos(mPanEnv());
            float s1 = (*this)() * mAmpEnv() * mAmp;
            float s2;
            mEnvFollow(s1);
            mPan(s1, s1, s2);
            io.out(0) += s1;
            io.out(1) += s2;
            // STFT for each notes
            if (stft(s1))
            { // Loop through all the frequency bins
                for (unsigned k = 0; k < stft.numBins(); ++k)
                {
                    // Here we simply scale the complex sample
                    spectrum[k] = tanh(pow(stft.bin(k).real(), 1.3));
                }
            }
        }
        if (mAmpEnv.done() && (mEnvFollow.value() < 0.001))
            free();
    }

    virtual void onTriggerOn() override
    {
        mAmpEnv.reset();
        timepose = 10;
        updateFromParameters();
        env.reset();
        delay.zero();
        mPanEnv.reset();
    }

    virtual void onTriggerOff() override
    {
        mAmpEnv.triggerRelease();
    }

    void updateFromParameters()
    {
        mPanEnv.levels(getInternalParameterValue("Pan1"),
                       getInternalParameterValue("Pan2"),
                       getInternalParameterValue("Pan1"));
        mPanRise = getInternalParameterValue("PanRise");
        delay.freq(getInternalParameterValue("frequency"));
        mAmp = getInternalParameterValue("amplitude");
        mAmpEnv.levels()[1] = 1.0;
        mAmpEnv.levels()[2] = getInternalParameterValue("sustain");
        mAmpEnv.lengths()[0] = getInternalParameterValue("attackTime");
        mAmpEnv.lengths()[3] = getInternalParameterValue("releaseTime");
        mPanEnv.lengths()[0] = mPanRise;
        mPanEnv.lengths()[1] = mPanRise;
    }
};

class MyApp : public App
{
public:
  // SynthGUIManager<SineEnv> synthManager {"synth8"};
  //    ParameterMIDI parameterMIDI;

  SynthGUIManager<DemoVoice> synthManager {"DemoVoice"};

  std::map<std::string, float> noteDictionary;
  std::map<int, std::string> numToNote;
  std::map<int, std::string> numToPentatonicRoot;
  std::map<std::string, vector<string>> pentatonicMajorMap;
  
  int measures = 0; //how many total measures, (4 beats per measure), measure lies in range [4,10]
  float measureTimeDuration = 0.0; //how long the measure lasts in seconds
  int beatsPerMeasure = 3; //how many beats per measure
  float timeIntervalBetweenBeats;
  // int bpm = getRandomInt(60,180); //randomly generate this on onInit() from [60, 180] //TODO
  int bpm;
  int pentatonicRootInt;
  std::string pentatonicRootString;

  Mesh mDemoVoice;
  vector<float> spectrum;
  
  int count = 0;

  float getFrequencyByOctave(std::string note, int octave) {
    return noteDictionary[note] * pow(2, octave);
  }

  int getRandomInt(int min, int max) {
    return rand() % (max - min + 1) + min;
  }

  float getAmplitudeFromTimePlayed(float time) { //a function that decreases the amplitude inversely with the time a note is played
    return 1/(time + 1);
  }

  virtual void onInit( ) override {
    imguiInit();
    navControl().active(false);  // Disable navigation via keyboard, since we
                              // will be using keyboard for note triggering
    // Set sampling rate for Gamma objects from app's audio
    gam::sampleRate(audioIO().framesPerSecond());

    // Initialize note dictionary with 0 octave frequencies
    noteDictionary["C"] = 16.351;
    noteDictionary["C# / Db"] = 17.324;
    noteDictionary["D"] = 18.354;
    noteDictionary["D# / Eb"] = 19.445;
    noteDictionary["E"] = 20.601;
    noteDictionary["F"] = 21.827;
    noteDictionary["F# / Gb"] =  23.124;
    noteDictionary["G"] = 24.499;
    noteDictionary["G# / Ab"] = 25.956;
    noteDictionary["A"] = 27.5;
    noteDictionary["A# / Bb"] = 29.135;
    noteDictionary["B"] = 30.868;

    // map numbers to notes
    numToNote[0] = "C";
    numToNote[1] = "C# / Db";
    numToNote[2] = "D";
    numToNote[3] = "D# / Eb";
    numToNote[4] = "E";
    numToNote[5] = "F";
    numToNote[6] = "F# / Gb";
    numToNote[7] = "G";
    numToNote[8] = "G# / Ab";
    numToNote[9] = "A";
    numToNote[10] = "A# / Bb";
    numToNote[11] = "B";

    //create all notes in a specific major pentatonic scale
    pentatonicMajorMap["C"] = {"C", "D", "E", "G", "A"};
    pentatonicMajorMap["D"] = {"D", "E", "F# / Gb", "A", "B"};
    pentatonicMajorMap["E"] = {"E", "F# / Gb", "G# / Ab", "B", "C# / Db"};
    pentatonicMajorMap["F"] = {"F", "G", "A", "C", "D"};
    pentatonicMajorMap["G"] = {"G", "A", "B", "D", "E"};
    pentatonicMajorMap["A"] = {"A", "B", "C# / Db", "E", "F# / Gb"};
    pentatonicMajorMap["G# / Ab"] = {"G# / Ab", "A# / Bb", "C", "D# / Eb", "F"};
    pentatonicMajorMap["B"] = {"B", "C# / Db", "D#", "F# / Gb", "G# / Ab"};

    //to map random number to a specific pentatonic root
    numToPentatonicRoot[0] = "C";
    numToPentatonicRoot[1] = "D";
    numToPentatonicRoot[2] = "E";
    numToPentatonicRoot[3] = "F";
    numToPentatonicRoot[4] = "G";
    numToPentatonicRoot[5] = "A";
    numToPentatonicRoot[6] = "B";

    //for piece timing and speed
    bpm = getRandomInt(60, 180);
    measures = getRandomInt(8, 25);
    measureTimeDuration = (60.0 / bpm) * (beatsPerMeasure);
    timeIntervalBetweenBeats = 1 / (bpm / 60.0);

    int pentatonicRootInt = getRandomInt(0, 6);
    pentatonicRootString = numToPentatonicRoot[pentatonicRootInt];

    cout << "Pentatonic Root: " << pentatonicRootString << endl;
    cout << "BPM: " << bpm << endl;
    cout << "MEASURES: " << measures << endl;
  }

    void onCreate() override {
        synthManager.synthRecorder().verbose(true);
    }

    void onSound(AudioIOData& io) override {
        synthManager.render(io);  // Render audio
    }

    void onAnimate(double dt) override {
        imguiBeginFrame();
        synthManager.drawSynthControlPanel();
        imguiEndFrame();
    }

    void onDraw(Graphics& g) override {
        g.clear();
        synthManager.render(g);

        // Draw GUI
        imguiDraw();
    }

    bool onKeyDown(Keyboard const& k) override {
        if (ParameterGUI::usingKeyboard()) {  // Ignore keys if GUI is using them
        return true;
        }
        if (k.shift()) {
        // If shift pressed then keyboard sets preset
        int presetNumber = asciiToIndex(k.key());
        synthManager.recallPreset(presetNumber);
        }
        if (k.key() == 'i') {
          startPiece();
        }
        if (k.key() == 'd') {
          createSnareDrumNote(100, 0.3, 0.1, 0.6, 0.0, 0, 0.2);
        }
        if (k.key() == 's') {
          createBassDrumNote(80, 0.3, 0.1, 0.6, 0.0, 0, 1.0);
        }
        if (k.key() == 'c') {
          createCymbalNote(1000, 0.3, 0.05, 2.0, 0.0, 0, 0.1);
        }
        if (k.key() == 'k') {
          int randomNote = getRandomInt(0, 4);
          createPluckedNote(noteDictionary[pentatonicMajorMap["G# / Ab"].at(randomNote)] * pow(2,4), 0.3, 0.001, 3.0, 0.0, 0, 0.1);
        }
        return true;
    }

    bool onKeyUp(Keyboard const& k) override {
        int midiNote = asciiToMIDI(k.key());
        if (midiNote > 0) {
        synthManager.triggerOff(midiNote);
        }
        return true;
    }

    void onExit() override { imguiShutdown(); }

    void createNote(float frequency, float amplitude, float attackTime, float releaseTime, float pan, double startTime, double duration) {
        auto* voice = synthManager.synth().getVoice<DemoVoice>();
        voice->setInternalParameterValue("frequency", frequency);
        voice->setInternalParameterValue("amplitude", amplitude);
        voice->setInternalParameterValue("attackTime", attackTime);
        voice->setInternalParameterValue("releaseTime", releaseTime);
        voice->setInternalParameterValue("pan", pan);
        synthManager.synthSequencer().addVoiceFromNow(voice, startTime, duration);
    }

    void createSnareDrumNote(float frequency, float amplitude, float attackTime, float releaseTime, float pan, double startTime, double duration) {
      auto* voice = synthManager.synth().getVoice<Snare>();
      voice->setInternalParameterValue("frequency", frequency);
      voice->setInternalParameterValue("amplitude", amplitude);
      voice->setInternalParameterValue("attackTime", attackTime);
      voice->setInternalParameterValue("releaseTime", releaseTime);
      voice->setInternalParameterValue("pan", pan);
      synthManager.synthSequencer().addVoiceFromNow(voice, startTime, duration);
    }

    void createBassDrumNote(float frequency, float amplitude, float attackTime, float releaseTime, float pan, double startTime, double duration) {
      auto* voice = synthManager.synth().getVoice<BassDrum>();
      voice->setInternalParameterValue("frequency", frequency);
      voice->setInternalParameterValue("amplitude", amplitude);
      voice->setInternalParameterValue("attackTime", attackTime);
      voice->setInternalParameterValue("releaseTime", releaseTime);
      voice->setInternalParameterValue("pan", pan);
      synthManager.synthSequencer().addVoiceFromNow(voice, startTime, duration);
    }

    void createCymbalNote(float frequency, float amplitude, float attackTime, float releaseTime, float pan, double startTime, double duration) {
      auto* voice = synthManager.synth().getVoice<Cymbal>();
      voice->setInternalParameterValue("frequency", frequency);
      voice->setInternalParameterValue("amplitude", amplitude);
      voice->setInternalParameterValue("attackTime", attackTime);
      voice->setInternalParameterValue("releaseTime", releaseTime);
      voice->setInternalParameterValue("pan", pan);
      synthManager.synthSequencer().addVoiceFromNow(voice, startTime, duration);
    }

    void createPluckedNote(float frequency, float amplitude, float attackTime, float releaseTime, float pan, double startTime, double duration) {
      auto* voice = synthManager.synth().getVoice<PluckedString>();
      voice->setInternalParameterValue("frequency", frequency);
      voice->setInternalParameterValue("amplitude", amplitude);
      voice->setInternalParameterValue("attackTime", attackTime);
      voice->setInternalParameterValue("releaseTime", releaseTime);
      voice->setInternalParameterValue("pan", pan);
      synthManager.synthSequencer().addVoiceFromNow(voice, startTime, duration);
    }

    void createAllDrumNotes() {
      int x = 0;
      float lastTime = 0.0;
      for(int i = 0; i < measures * beatsPerMeasure; i++) {
        createSnareDrumNote(100, 0.3, 0.1, 0.6, 0.0,i * timeIntervalBetweenBeats, 0.3);   
        // createBassDrumNote(80, 0.3, 0.1, 0.6, 0.0, i * timeIntervalBetweenBeats, 0.3);     
      }
    }
    void createAllMelodyNotes() {
      int genOrNo = getRandomInt(0, 1);
      for(int i = 0; i < measures * beatsPerMeasure; i++) {
        if(genOrNo == 1) {
          int randomNote = getRandomInt(0, 4);
          createNote(noteDictionary[pentatonicMajorMap[pentatonicRootString].at(randomNote)] * pow(2,4), 0.3, 0.1, 0.6, 0.0, i * timeIntervalBetweenBeats, 0.3);
        } 
        genOrNo = getRandomInt(0, 1);
      }
    }
    void createBackgroundChords() {
      for(int i = 0; i < measures; i++) {
        createNote(noteDictionary[pentatonicMajorMap[pentatonicRootString].at(i % 5)] * pow(2, 3), 0.1, 0.1, 0.1, 0.0, i * timeIntervalBetweenBeats * beatsPerMeasure, 3);
      }
    }

    void startPiece() {
        createAllDrumNotes();
        createAllMelodyNotes();
        createBackgroundChords();
    }
};

int main() {
  MyApp app;

  // Set up audio
  app.configureAudio(48000., 512, 2, 0);

  app.start();
}