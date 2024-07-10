#include <bits/stdc++.h>
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

const double PI = 3.14159265358979323846;

// Complex number struct
struct Complex {
    double real;
    double imag;

    Complex(double r = 0, double i = 0) : real(r), imag(i) {}

    Complex operator+(const Complex& other) const {
        return Complex(real + other.real, imag + other.imag);
    }

    Complex operator-(const Complex& other) const {
        return Complex(real - other.real, imag - other.imag);
    }

    Complex operator*(const Complex& other) const {
        return Complex(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }

    Complex operator*(double scalar) const {
        return Complex(real * scalar, imag * scalar);
    }

    Complex conj() const {
        return Complex(real, -imag);
    }
};

// Function to read a .wav file and return signal data and sample rate
bool readWavFile(const char* filename, vector<double>& signal, int& sampleRate) {
    ifstream file(filename, ios::binary);
    if (!file) {
        cerr << "Error opening file " << filename << endl;
        return false;
    }

    // Read header (assuming 16-bit PCM WAV format)
    char header[44];
    file.read(header, 44);

    // Extract sample rate from header
    sampleRate = *reinterpret_cast<int*>(header + 24);
    
    


    // Read data
    short buffer;
    while (file.read(reinterpret_cast<char*>(&buffer), sizeof(buffer))) {
        signal.push_back(static_cast<double>(buffer) / 32768.0); // Convert to range [-1, 1]
        cout<<static_cast<double>(buffer) / 32768.0<<endl;
    }

    return true;
}

// Function to write a .wav file from signal data
bool writeWavFile(const char* filename, const vector<double>& signal, int sampleRate) {
    ofstream file(filename, ios::binary);
    if (!file) {
        cerr << "Error opening file " << filename << endl;
        return false;
    }

    // Prepare header for 16-bit PCM WAV file
    const int numSamples = signal.size();
    const int byteRate = sampleRate * sizeof(short);
    const int blockAlign = sizeof(short);

    file << "RIFF";
    file.write(reinterpret_cast<const char*>(&numSamples), 4); // Chunk size
    file << "WAVEfmt ";
    int subchunk1Size = 16;
    file.write(reinterpret_cast<const char*>(&subchunk1Size), 4); // Subchunk1Size
    short audioFormat = 1; // PCM
    file.write(reinterpret_cast<const char*>(&audioFormat), 2);
    short numChannels = 1; // Mono
    file.write(reinterpret_cast<const char*>(&numChannels), 2);
    file.write(reinterpret_cast<const char*>(&sampleRate), 4);
    file.write(reinterpret_cast<const char*>(&byteRate), 4);
    file.write(reinterpret_cast<const char*>(&blockAlign), 2);
    short bitsPerSample = 16;
    file.write(reinterpret_cast<const char*>(&bitsPerSample), 2);
    file << "data";
    int subchunk2Size = numSamples * sizeof(short);
    file.write(reinterpret_cast<const char*>(&subchunk2Size), 4);

    // Write sample data
    for (double sample : signal) {
        short intSample = static_cast<short>(sample * 32767);
        file.write(reinterpret_cast<const char*>(&intSample), sizeof(intSample));
    }

    return true;
}

// Perform the FFT algorithm
void fft(vector<Complex>& a) {
    int n = a.size();
    if (n <= 1) return;

    // Divide
    vector<Complex> a0(n / 2), a1(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        a0[i] = a[i * 2];
        a1[i] = a[i * 2 + 1];
    }
    fft(a0);
    fft(a1);

    // Conquer
    double angle = 2 * PI / n;
    Complex w(1), wn(cos(angle), sin(angle));
    for (int i = 0; i < n / 2; ++i) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        if (i + n / 2 < n) {
            w = w * wn;
        }
    }
}

// Perform the inverse FFT algorithm
void ifft(vector<Complex>& a) {
    int n = a.size();
    if (n <= 1) return;

    // Conjugate the complex numbers
    for (Complex& x : a) {
        x = x.conj();
    }

    // Forward FFT
    fft(a);

    // Conjugate the complex numbers again
    for (Complex& x : a) {
        x = x.conj();
    }

    // Divide by the number of samples
    for (Complex& x : a) {
        x = x * (1.0 / n);
    }
}

// Apply a low-pass filter by zeroing out high frequencies
void lowPassFilter(vector<Complex>& fftData, int sampleRate, double cutoff) {
    int n = fftData.size();
    double nyquist = sampleRate / 2.0;

    for (int i = 0; i < n; ++i) {
        double freq = i * nyquist / (n / 2);
        if (freq > cutoff) {
            fftData[i] = Complex(0, 0);
        }
    }
}

int main() {
    const char* inputFile = "C:\\Users\\abhis\\Downloads\\sweep.wav";
    const char* outputFile = "C:\\Users\\abhis\\Downloads\\filtered_output11.wav";

    vector<double> signal;
    int sampleRate;

    // Read the .wav file
    if (!readWavFile(inputFile, signal, sampleRate)) {
        return 1;
    }

    // Perform FFT
    vector<Complex> fftData(signal.begin(), signal.end());
    fft(fftData);

    // Apply a low-pass filter
    double cutoff = 1000.0; // Cutoff frequency in Hz
    lowPassFilter(fftData, sampleRate, cutoff);

    // Perform inverse FFT
    ifft(fftData);

    // Extract the real part of the filtered signal
    vector<double> filteredSignal(signal.size());
    for (size_t i = 0; i < filteredSignal.size(); ++i) {
        filteredSignal[i] = fftData[i].real;
    }

    //Write the filtered signal to a new .wav file
    if (!writeWavFile(outputFile, filteredSignal, sampleRate)) {
        return 1;
    }

    cout << "Filtered signal saved to " << outputFile << endl;

    return 0;
}
