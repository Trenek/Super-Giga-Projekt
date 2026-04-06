#include <unistd.h>

#include "draw.hpp"

static void setGNUPlot(int id, struct thing thing) {
    fprintf(thing.fileDesc, "set term qt %d size 800,600\n", id);
    fprintf(thing.fileDesc, "set title '%s'\n", thing.name);
    fprintf(thing.fileDesc, "set xlabel '%s'\n", thing.xName);
    fprintf(thing.fileDesc, "set ylabel '%s'\n", thing.yName);
    fprintf(thing.fileDesc, "set grid\n");

    fprintf(thing.fileDesc, "pause 1\n");
    fprintf(thing.fileDesc, "plot \"%s\" with lines\n", thing.file);
    fprintf(thing.fileDesc, "while (1) {\n");
    fprintf(thing.fileDesc, "    pause 1\n");
    fprintf(thing.fileDesc, "    replot\n");
    fprintf(thing.fileDesc, "}\n");
    fflush(thing.fileDesc);
}

static void cleanup(std::vector<struct thing> thing) {
    for (size_t i = 0; i < thing.size(); i += 1) {
        std::remove(thing[i].file);
    }
    sleep(4);
}

static void createPlot(std::vector<struct thing> thing) {
    for (size_t i = 0; i < thing.size(); i += 1) {
        thing[i].fileDesc = popen("gnuplot", "w");

        setGNUPlot(i, thing[i]);
    }
};

static void destroyPlot(std::vector<struct thing> thing) {
    for (size_t i = 0; i < thing.size(); i += 1) {
        fclose(thing[i].fileDesc);
    }
}

gnuPlotManager::gnuPlotManager(std::vector<struct thing> &&array) : array(std::move(array)) {
    cleanup(this->array);

    for (size_t i = 0; i < this->array.size(); i += 1) {
        this->array[i].fileDesc = fopen(this->array[i].file, "w");
    }

    createPlot(this->array);
}

gnuPlotManager::~gnuPlotManager() {
    destroyPlot(this->array);
}
