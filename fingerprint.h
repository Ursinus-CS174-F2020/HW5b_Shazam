#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include "hashable.h"
#include "hashtable.h"
#include "BMP.h"
#include "audio.h"
#include "dsp.h"
#include "dspviz.h"
#include <vector>
#include <fstream>
#include <string.h>
#include <string>
#include <sstream>
using namespace std;


class Anchor {
    public:
        int freq; // Frequency of occurrence
        int win; // Window number in the audio
        Anchor(int freq, int win);
        Anchor() {}

        /**
         * Find all anchors that are within the target zone of this anchor
         * @param others A list of anchors to consider
         * @param dcenter Offset in time of center of window
         * @param width Half-width of window in time
         * @param height Half-height of window in frequency
         */
        vector<Anchor>* searchTargetZone(vector<Anchor>* others, int dcenter, int width, int height);
};

class Fingerprint: public Hashable {
    public:
        Anchor a1; // First local max
        Anchor a2; // Second local max (assumed to be after the first)
        string name; // Song name associated to this fingerprint
        Fingerprint(){}
        void setAnchors(Anchor* a1, Anchor* a2);
        Fingerprint(Anchor* a1, Anchor* a2);
        size_t getHash();
        bool equals(Hashable* other);
        Cloneable* clone();
};


/**
 * Figure out whether a particular point is a max
 * @param S Spectrogram
 * @param i Frequency index
 * @param j Time index
 * @param maxFreq Maximum frequency to consider
 * @param nwin Number of windows in time
 * @param freqWin Half-length of window in frequency
 * @param timeWin Half-length of window in time
 * @return true if this element is greater than all elements in the window, false otherwise
 */
bool isMax(double** S, int i, int j, int maxFreq, int nwin, int freqWin, int timeWin);


/**
 * Figure out whether a particular point is a max
 * @param S Spectrogram
 * @param maxFreq Maximum frequency to consider
 * @param nwin Number of windows in time
 * @param freqWin Half-length of window in frequency
 * @param timeWin Half-length of window in time
 * @param thresh Minimum value in spectrogram to consider
 * @return true if this element is greater than all elements in the window, false otherwise
 */
vector<Anchor>* findAnchors(double** S, int maxFreq, int nwin, int freqWin, int timeWin, double thresh);


/**
 * Compute the fingerprints in a set of anchors
 * @param anchors A list of anchors
 * @param dcenter Offset in time of center of window
 * @param width Half-width of window in time
 * @param height Half-height of window in frequency
 */
vector<Fingerprint>* getFingerprints(vector<Anchor>* anchors, int dcenter, int width, int height);


/**
 * Plot the anchors as dots on a canvas
 * @param canvas BMP canvas to which to draw dots
 * @param anchors List of anchors to plot
 * @param scale How much the canvas is scaled
 */
void plotAnchors(BMP& canvas, vector<Anchor>* anchors, int scale);

/**
 * Plot the fingerprints as lines on a canvas
 * @param canvas BMP canvas to which to draw dots
 * @param fingerprints List of fingerprints to plot
 * @param scale How much the canvas is scaled
 */
void plotFingerprints(BMP& canvas, vector<Fingerprint>* fingerprints, int scale);

typedef struct SongInfo {
    string title;
    string artist;
    void printInfo() {
        cout << artist << ": " << title << endl;
    }
} songInfo;

class SongData {
    public:
        songInfo info;
        vector<Fingerprint>* fingerprints;
        bool outputImage;
        const char* impath;
        bool outputFingerprints;
        const char* fpath;
        bool audioPathSpecified;
        const char* audioPath;
        short* audio;
        size_t N;
        ~SongData();
        SongData(){}
        /**
         * Update the audio path of this song
         * @param audioPath String path to audio
         */
        void setAudioPath(const char* audioPath);
        /**
         * Parse command line arguments for audio path and other options like 
         * where to output bitmap showing fingerprints
         * @param argc Number of arguments
         * @param argv Command line arguments
         */
        void parseArgs(int argc, char** argv);
        /**
         * Construct a song from command line arguments
         * @param argc Number of arguments
         * @param argv Command line arguments
         */
        SongData(int argc, char** argv);
        /**
         * Load in audio and construct fingerprints.  This method
         * assumes that a path to the audio exists.  It also plots
         * the fingerprints and saves them to a file if those options
         * were specified
         */
        void processAudio();

        /**
         * Query the database for the best matching song
         * @param db Hash table with all of the fingerprints
         * @param songs Array of info for songs in the database
         * @param NSongs Number of songs in the database
         */
        songInfo queryDatabase(HashTable* db, songInfo* songs, int NSongs);
};

/**
 * A struct to store information to associate a fingerprint in
 * the database to its song
 */
typedef struct FingerprintInfo {
    int hashCode;
    int offset;
    int index; // Index in database
    string artist;
    string title;
} fingerprintInfo;

/**
 * A class to wrap around a vector for storing fingerprints of the same
 * type for different songs in a list.  The wrapper is cloneable so
 * this list can be used as a value in the HashMap
 */
class FingerprintList:public Cloneable {
    public:
        vector<fingerprintInfo> list;
        Cloneable* clone() {
            FingerprintList* L = new FingerprintList();
            for (size_t i = 0; i < list.size(); i++) {
                L->list.push_back(list.at(i));
            }
            return L;
        }
        size_t size() {
            return list.size();
        }
        fingerprintInfo at(size_t i) {
            return list.at(i);
        }
        void push_back(fingerprintInfo fingerprint) {
            list.push_back(fingerprint);
        }
};

/**
 * Load in a database of fingerprints of different songs
 * @param N Number of songs to load in
 * @param NBins Number of bins in the hash table
 * @param songs Pointer to the array that will hold the songs
 */
HashTable* loadDatabase(int N, int NBins, songInfo** songs);


#endif