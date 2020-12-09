#include "fingerprint.h"

Anchor::Anchor(int freq, int win) {
    this->freq = freq;
    this->win = win;
}

/**
 * Find all anchors that are within the target zone of this anchor
 * @param others A list of anchors to consider
 * @param dcenter Offset in time of center of window
 * @param width Half-width of window in time
 * @param height Half-height of window in frequency
 */
vector<Anchor>* Anchor::searchTargetZone(vector<Anchor>* others, int dcenter, int width, int height) {
    vector<Anchor>* matches = new vector<Anchor>();
    // TODO: Fill this in
    return matches;
}

void Fingerprint::setAnchors(Anchor* a1, Anchor* a2) {
    this->a1.freq = a1->freq;
    this->a1.win = a1->win;
    this->a2.freq = a2->freq;
    this->a2.win = a2->win;
}
Fingerprint::Fingerprint(Anchor* a1, Anchor* a2) {
    setAnchors(a1, a2);
}
size_t Fingerprint::getHash() {
    // TODO: Fill this in
    return 0; // This is a dummy value
}
bool Fingerprint::equals(Hashable* other) {
    Fingerprint* otherp = (Fingerprint*)other;
    return getHash() == otherp->getHash();
}
Cloneable* Fingerprint::clone() {
    return new Fingerprint(&a1, &a2);
}

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
bool isMax(double** S, int i, int j, int maxFreq, int nwin, int freqWin, int timeWin) {
    bool max = true;
    // TODO: Fill this in
    return max;
}

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
vector<Anchor>* findAnchors(double** S, int maxFreq, int nwin, int freqWin, int timeWin, double thresh) {
    vector<Anchor>* anchors = new vector<Anchor>();
    // TODO: Fill this in
    return anchors;
}

/**
 * Compute the fingerprints in a set of anchors
 * @param anchors A list of anchors
 * @param dcenter Offset in time of center of window
 * @param width Half-width of window in time
 * @param height Half-height of window in frequency
 */
vector<Fingerprint>* getFingerprints(vector<Anchor>* anchors, int dcenter, int width, int height) {
    vector<Fingerprint>* fingerprints = new vector<Fingerprint>();
    // TODO: Fill this in
    return fingerprints;
}

/**
 * Plot the anchors as dots on a canvas
 * @param canvas BMP canvas to which to draw dots
 * @param anchors List of anchors to plot
 * @param scale How much the canvas is scaled
 */
void plotAnchors(BMP& canvas, vector<Anchor>* anchors, int scale) {
    for (size_t i = 0; i < anchors->size(); i++) {
        Anchor a = anchors->at(i);
        canvas.plotCircle(a.win*scale+scale/2, a.freq*scale+scale/2, scale, 0xFF, 0xFF, 0);
    }
}

/**
 * Plot the fingerprints as lines on a canvas
 * @param canvas BMP canvas to which to draw dots
 * @param fingerprints List of fingerprints to plot
 * @param scale How much the canvas is scaled
 */
void plotFingerprints(BMP& canvas, vector<Fingerprint>* fingerprints, int scale) {
    for (size_t i = 0; i < fingerprints->size(); i++) {
        Anchor* a1 = &fingerprints->at(i).a1;
        Anchor* a2 = &fingerprints->at(i).a2;
        int x1 = a1->win*scale+scale/2;
        int x2 = a2->win*scale+scale/2;
        int y1 = a1->freq*scale+scale/2;
        int y2 = a2->freq*scale+scale/2;
        canvas.plotLine(x1, y1, x2, y2, 1, 0, 0, 0xFF);
    }
}



SongData::~SongData() {
    if (fingerprints != NULL) {
        delete fingerprints;
    }
}

/**
 * Update the audio path of this song
 * @param audioPath String path to audio
 */
void SongData::setAudioPath(const char* audioPath) {
    this->audioPath = audioPath;
    audioPathSpecified = true;
}

/**
 * Parse command line arguments for audio path and other options like 
 * where to output bitmap showing fingerprints
 * @param argc Number of arguments
 * @param argv Command line arguments
 */
void SongData::parseArgs(int argc, char** argv) {
    outputImage = false;
    outputFingerprints = false;
    audioPathSpecified = false;
    argv++, argc--;
    info.title = string("untitled");
    info.artist = string("anon");
    while (argc > 0) {
        if ((*argv)[0] == '-') {
            if (strcmp(*argv, "-audio") == 0) {
                argv++; argc--;
                setAudioPath((const char*)*argv);
            }
            else if (!strcmp(*argv, "-title")) { 
                argv++; argc--;
                info.title = string((const char*)*argv);
            }
            else if (!strcmp(*argv, "-artist")) {
                argv++; argc--;
                info.artist = string((const char*)*argv);
            }
            else if (!strcmp(*argv, "-fpath")) {
                argv++; argc--;
                fpath = (const char*)*argv;
                outputFingerprints = true;
            }
            else if (!strcmp(*argv, "-impath")) {
                argv++; argc--;
                impath = (const char*)*argv;
                outputImage = true;
            }
            else { 
                fprintf(stderr, "Invalid option: %s\n", *argv);
            }
        }
        argv++, argc--; 
    }
}

/**
 * Construct a song from command line arguments to an audio file path
 * and other options like where to output bitmap showing fingerprints
 * @param argc Number of arguments
 * @param argv Command line arguments
 */
SongData::SongData(int argc, char** argv) {
    parseArgs(argc, argv);
    fingerprints = NULL;
}


/**
 * Load in audio and construct fingerprints.  This method
 * assumes that a path to the audio exists.  It also plots
 * the fingerprints and saves them to a file if those options
 * were specified
 */
void SongData::processAudio() {
    if (!audioPathSpecified) {
        printf("Error: Need to at least specify audio path\n");
        return;
    }
    // Assuming audio is sampled at 22050hz
    short* audio = getAudio((const char*)audioPath, &N);
    int win = 2048;
    int maxFreq = 128;
    int hop = 1024;
    DSP dsp(win);
    int nwin = 0;
    printf("Getting spectrogram...\n");
    double** specgram = dsp.specgram(audio, N, win, hop, true, &nwin);
    printf("Finished getting spectrogram\n");
    vector<Anchor>* anchors = findAnchors(specgram, maxFreq, nwin, 8, 3, 0);
    fingerprints = getFingerprints(anchors, 86, 50, 21);
    if (outputFingerprints) {
        ofstream fout;
        fout.open(fpath);
        fout << info.artist << "\n" << info.title;
        for (size_t i = 0; i < fingerprints->size(); i++) {
            fout << "\n" << fingerprints->at(i).getHash();
            fout << " " << fingerprints->at(i).a1.win;
        }
        fout.close();
    }
    if (outputImage) {
        BMP canvas = plotSpectrogram(specgram, maxFreq, nwin, 4, 4);
        plotFingerprints(canvas, fingerprints, 4);
        plotAnchors(canvas, anchors, 4);
        canvas.write(impath);
    }
    deleteSpecgram(specgram, win);
    delete[] audio;
    delete anchors;
}

/**
 * Query the database for the best matching song
 * @param db Hash table with all of the fingerprints
 * @param songs Array of info for songs in the database
 * @param NSongs Number of songs in the database
 */
songInfo SongData::queryDatabase(HashTable* h, songInfo* songs, int NSongs) {
    vector<HashTable*> histograms;
    int maxCount = 0;
    int maxIndex = 0;
    for (int i = 0; i < NSongs; i++) {
        HashTable* h = new HashTable(100);
        histograms.push_back(h);
    }

    
    // TODO: Fill this in
    
    
    
    for (int i = 0; i < NSongs; i++) {
        delete histograms.at(i);
    }
    return songs[maxIndex];
}


/**
 * Load in a database of fingerprints of different songs
 * @param N Number of songs to load in
 * @param NBins Number of bins in the hash table
 * @param songs Pointer to the array that will hold the songs
 */
HashTable* loadDatabase(int N, int NBins, songInfo** songs) {
    HashTable* h = new HashTable(NBins);
    *songs = new songInfo[N];
    for (int i = 0; i < N; i++) {
        stringstream ss;
        ss << "Database/" << i << ".txt";
        ifstream file(ss.str());
        string artist;
        string title;
        int lineIdx = 0;
        if (file.is_open()) {
            string line;
            while(getline(file, line)) {
                if (lineIdx == 0) {
                    artist = line;
                }
                else if (lineIdx == 1) {
                    title = line;
                    cout << "Loading " << artist << ": " << title << "\n";
                    (*songs)[i].artist = artist;
                    (*songs)[i].title = title;
                }
                else {
                    int pos = line.find(" ");
                    int f = stoi(line.substr(0, pos));
                    line.erase(0, pos+1);
                    int offset = stoi(line);
                    Fingerprint fingerprint;
                    fingerprint.a1.win = offset;
                    fingerprint.a1.freq = f & 0xFF;
                    fingerprint.a2.win = offset + ((f >> 16) & 0xFF);
                    fingerprint.a2.freq = (f >> 8) & 0xFF;
                    fingerprintInfo info;
                    info.hashCode = f;
                    info.offset = offset;
                    info.artist = artist;
                    info.title = title;
                    info.index = i;
                    if (!h->containsKey(&fingerprint)) {
                        FingerprintList temp;
                        h->put(&fingerprint, &temp);
                    }
                    FingerprintList* flist = (FingerprintList*)h->get(&fingerprint);
                    flist->push_back(info);
                }
                lineIdx++;
            }
        }
        else {
            cout << "Warning: Could not open file " << ss.str() << "\n";
        }
    }
    return h;
}
