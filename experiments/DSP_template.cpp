// Joel A. Jaffe 2024-04-01
// Basic 2D AlloApp with Oscilliscope

#include "al/app/al_App.hpp"
#include "al/graphics/al_Mesh.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/math/al_Random.hpp"
using namespace al;

#include <iostream>
using namespace std;

// Useful functions for audio
float dBtoA (float dBVal) {return powf(10.f, dBVal / 20.f);}
float ampTodB (float ampVal) {return 20.f * log10f(fabs(ampVal));}
float mToF (int midiVal) {return 440.f * powf(2.f, (midiVal - 69) / 12.f);}
int fToM (float freq) {return 12.f * log2f(freq / 440.f) + 69;}

// Oscilliscope class that inherits from Mesh
class Oscilliscope : public Mesh {
public:
  Oscilliscope (int samplerate) : bufferSize(samplerate) { 
    this->primitive(Mesh::LINE_STRIP);
    for (int i = 0; i < bufferSize; i++) {
      this->vertex((i / static_cast<float>(bufferSize)) * 2.f - 1.f, 0);
      buffer.push_back(0.f);
    }
  }

  void writeSample (float sample) { // call inside while(io()) loop inside onAudio() block
    for (int i = 0; i < bufferSize - 1; i++) {
      buffer[i] = buffer[i + 1];
    }
    buffer[bufferSize - 1] = sample;
  }

  void update() { // call inside onAnimate() block
    for (int i = 0; i < bufferSize; i++) {
      this->vertices()[i][1] = buffer[i];
    }
  }
    
protected:
  int bufferSize;
  vector<float> buffer;
};

// Barebones app that displays an oscilliscope 
struct DSP_Template : public App {
  Parameter volControl{"volControl", "", 0.f, -96.f, 6.f};
  Parameter rmsMeter{"rmsMeter", "", -96.f, -96.f, 0.f};
  ParameterBool audioOutput{"audioOutput", "", false, 0.f, 1.f};
  Oscilliscope myScope{static_cast<int>(AudioIO().framesPerSecond())};

  void onInit() {
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(volControl); // add parameter to GUI
    gui.add(rmsMeter);
    gui.add(audioOutput); 
  }

  void onCreate() {}

  void onAnimate(double dt) {
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
      io.out(0) = rnd::uniformS() * volFactor * audioOutput; // write white-noise to L channel
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
    g.clear(0); // Draws background black
    g.color(1); // Draws oscilliscope white
    g.camera(Viewpoint::IDENTITY); // Locks camera 
    g.draw(myScope); // Draws our mesh
  }
};
  
int main() {
  DSP_Template app;

  // Allows for manual declaration of input and output devices, 
  // but causes unpredictable behavior. Needs investigation.
  app.audioIO().deviceIn(AudioDevice("MacBook Pro Microphone")); // Configure your Audio I
  app.audioIO().deviceOut(AudioDevice("MacBook Pro Speakers")); //            and O devices
  app.configureAudio(44100, 128, app.audioIO().channelsOutDevice(), app.audioIO().channelsInDevice());

  /* 
  // Declaration of AudioDevice using aggregate device
  AudioDevice alloAudio = AudioDevice("AlloAudio");
  alloAudio.print();
  app.player.rate(1.0 / alloAudio.channelsOutMax());
  app.configureAudio(alloAudio, 44100, 128, alloAudio.channelsOutMax(), 2);
  */

  app.start();
}