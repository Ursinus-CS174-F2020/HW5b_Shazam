#include "stdio.h"
#include "audio.h"
#include <vector>
#include <string>
#include <fstream>
#include <cstdint>

using namespace std;


/**
 * Get the header for a wave file
 * @param filename Path to file
 * @return Wave header
 */
wav_hdr getAudioHeader(char* filename) {
    wav_hdr wavHeader;
    FILE* wavFile = fopen(filename, "r");
    if (wavFile == NULL) {
        fprintf(stderr, "Unable to open wave file: %s\n", filename);
    }
    size_t bytesRead = fread(&wavHeader, 1, sizeof(wav_hdr), wavFile);
    if (bytesRead != sizeof(wav_hdr)) {
        fprintf(stderr, "Unable to read audio header");
    }
    return wavHeader;
}

/**
 * Get the samples of an audio file
 * @param filename Path to audio file
 * @param N Pointer to file size
 * @return Short array of audio samples
 */
float* getAudio(char* filename, size_t* N) {
    short* ret = NULL;
    wav_hdr wavHeader;
    FILE* wavFile = fopen(filename, "r");
    if (wavFile == NULL) {
        fprintf(stderr, "Unable to open wave file: %s\n", filename);
    }
    else {
        size_t bytesRead = fread(&wavHeader, 1, sizeof(wav_hdr), wavFile);
        if (bytesRead != sizeof(wav_hdr)) {
            fprintf(stderr, "Unable to read audio header");
        }
        else if (wavHeader.NumOfChan != 1) {
            fprintf(stderr, "This library only handles mono audio");
        }
        else if (wavHeader.bitsPerSample != 16) {
            fprintf(stderr, "This library only handles 16-bit audio, but audio is %i-bit", wavHeader.bitsPerSample);
        }
        else {
            vector<float> samples;
            short chunk[AUDIO_CHUNK_SIZE];
            size_t bytesRead;
            *N = 0;
            do {
                bytesRead = fread(chunk, 1, AUDIO_CHUNK_SIZE, wavFile);
                for (size_t i = 0; i < bytesRead/2; i++) {
                    samples.push_back(chunk[i]);
                }
                *N += bytesRead/2;
            }
            while (bytesRead > 0);
            ret = new short[*N];
            for (size_t i = 0; i < *N; i++) {
                ret[i] = samples.at(i);
            }
        }

    }
    fclose(wavFile);
    return ret;
}

/**
 * Write 16-bit audio to a file
 * @param filename Path to output file
 * @param arr Audio samples to write (of type short)
 * @param N Number of samples
 * @param waveHeader Header to wave file
 */
void writeAudio(char* filename, short* arr, size_t N, wav_hdr wavHeader) {
    FILE* wavFile = fopen(filename, "w");
    fwrite(&wavHeader, 1, sizeof(wav_hdr), wavFile);
    fwrite(arr, 1, N*2, wavFile);
    fclose(wavFile);
}