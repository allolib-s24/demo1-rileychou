// Joel A. Jaffe 2024-04-10
// Basic 2D AlloApp with Oscilliscope, Phase Accumulator and Sine Osc 

#include "al/app/al_App.hpp"
#include "al/graphics/al_Mesh.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/math/al_Random.hpp"
using namespace al;

#include <iostream>
using namespace std;

float dBtoA (float dBVal) {return powf(10.f, dBVal / 20.f);}
float ampTodB (float ampVal) {return 20.f * log10f(fabs(ampVal));}
float mToF (int midiVal) {return 440.f * powf(2.f, (midiVal - 69) / 12.f);}
int fToM (float freq) {return 12.f * log2f(freq / 440.f) + 69;}

// Oscilliscope that inherits from mesh 
class Oscilliscope : public Mesh {
public:
  Oscilliscope (int samplerate) : bufferSize(samplerate) {
    this->primitive(Mesh::LINE_STRIP);
    for (int i = 0; i < bufferSize; i++) {
      this->vertex((i / static_cast<float>(bufferSize)) * 2.f - 1.f, 0);
      buffer.push_back(0.f);
    }
  }

  void writeSample (float sample) {
    for (int i = 0; i < bufferSize - 1; i++) {
      buffer[i] = buffer[i + 1];
    }
    buffer[bufferSize - 1] = sample;
  }

  void update() {
    for (int i = 0; i < bufferSize; i++) {
      this->vertices()[i][1] = buffer[i];
    }
  }
    
protected:
  int bufferSize;
  vector<float> buffer;
};

// Basic Phase Accumulator
// Generates unipolar ramp wave from 0 to 1 with adjustable frequency
class Phasor {
public:
  Phasor (int samprate) : sampleRate(samprate) {}

  virtual void setSampleRate (int samprate) {
    sampleRate = samprate;
    phaseIncrement = frequency / static_cast<float> (sampleRate);
  }

  virtual void setFrequency (float freq) {
    frequency = freq;
    phaseIncrement = frequency / static_cast<float> (sampleRate);
  }

  virtual float processSample() {
    float phaseIncrement = frequency / static_cast<float>(sampleRate);
    phase += phaseIncrement;
    phase = fmod(phase, 1.f);
    return phase;
  }

  virtual float setPhase (float ph) {phase = ph;}
  virtual float getPhase () {return phase;}

protected:
  float phase = 0.f;
  int sampleRate = 44100;
  float frequency = 1.f;
  float phaseIncrement = frequency / static_cast<float>(sampleRate);
};

// Sine Oscillator using Nth-order Taylor Series approximation, N adjustable
// as seen in Gamma https://github.com/LancePutnam/Gamma 
class SinOsc : public Phasor {
public: 
  SinOsc (int samprate) : Phasor(samprate) {}

  float processSample() override {
    phase += phaseIncrement;
    phase = fmod(phase, 1.f);
   return taylorNSin(phase * -twoPi + pi, N); // map phase to range (-pi,pi), calculate sin
  }

  float setOrder (int order) {N = order;}

protected:
  int N = 11; // 11 is good
  const float pi = static_cast<float>(M_PI);
  const float twoPi = 2.f * pi;

  int factorial (int x) {
    int output = 1;
    int n = x;
    while (n > 1) {
      output *= n;
      n -= 1;
    }
    return output;
  }

  float taylorNSin (float x, int order) {
    float output = x;
    int n = 3;
    int sign = -1;
    while (n <= order) {
     output += sign * (powf(x, n) / factorial(n));
      sign *= -1;
      n += 2;
    }
    return output;
  }
};

// app struct
struct DSP_Template : public App {
  Parameter volControl{"volControl", "", 0.f, -96.f, 6.f};
  Parameter rmsMeter{"rmsMeter", "", -96.f, -96.f, 0.f};
  Parameter myFreq{"myFreq", "", 220.f, 1.f, 10000.f};
  ParameterBool audioOutput{"audioOutput", "", false, 0.f, 1.f};
  Oscilliscope myScope{static_cast<int>(AudioIO().framesPerSecond())};
  SinOsc myOsc{static_cast<int>(AudioIO().framesPerSecond())};

  void onInit() {
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(volControl); // add parameter to GUI
    gui.add(rmsMeter);
    gui.add(audioOutput); 
    gui.add(myFreq);
  }

  void onCreate() {
    myOsc.setFrequency(220.f);
  }

  float last = 0.f;

  void onAnimate(double dt) {

    float changedFreq = (myFreq + last) / 2;
    last = changedFreq;
    myScope.update();
  }

  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == 'm') { // <- on m, muteToggle
      audioOutput = !audioOutput;
      cout << "Mute Status: " << audioOutput << endl;
    }
    return true;
  }

  void onSound(AudioIOData& io) override {
    // audio throughput and analysis
    float bufferPower = 0;
    float volFactor = dBtoA(volControl);
    while(io()) { 
      io.out(0) = myOsc.processSample() * volFactor * audioOutput; // write white-noise to L channel
      io.out(1) = io.out(0); // copy L channel to R channel
      myScope.writeSample((io.out(0) + io.out(1)) / 2.f); // write samples to osc
    }

    // feed to analysis buffer
    for (int channel = 0; channel < io.channelsIn(); channel++){
        bufferPower += powf(io.out(channel), 2);
    }
    bufferPower /= io.framesPerBuffer();
    rmsMeter = ampTodB(bufferPower);

    // overload detector
      if (io.out(0) > 1.f || io.out(1) > 1.f) {
        cout << "CLIP!" << endl;
    }
  }

  void onDraw(Graphics &g) {
    g.clear(0);
    g.color(1);
    g.camera(Viewpoint::IDENTITY); // Ortho [-1:1] x [-1:1]
    g.draw(myScope);
  }
};
  
int main() {
  DSP_Template app; // instance of our app

  // Allows for manual declaration of input and output devices, 
  // but causes unpredictable behavior. Needs investigation.
  app.audioIO().deviceIn(AudioDevice("MacBook Pro Microphone"));
  app.audioIO().deviceOut(AudioDevice("MacBook Pro Speakers"));
  app.configureAudio(44100, 128, app.audioIO().channelsOutDevice(), app.audioIO().channelsInDevice());

  /* 
  // Declaration of AudioDevice using aggregate device
  AudioDevice alloAudio = AudioDevice("AlloAudio");
  alloAudio.print();
  app.player.rate(1.0 / alloAudio.channelsOutMax());
  app.configureAudio(alloAudio, 44100, 128, alloAudio.channelsOutMax(), 2);
  */

  app.start();
  return 0;
}