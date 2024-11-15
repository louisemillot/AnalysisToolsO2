#ifndef GLOBAL_SETTINGS_H
#define GLOBAL_SETTINGS_H

// Bin edge control
const double GLOBAL_epsilon = 0.0000000000000000001;
const std::array<std::array<float, 2>, 2> drawnWindowAuto = {{{-999, -999}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}
const std::array<std::array<float, 2>, 2> legendPlacementAuto = {{{-999, -999}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}
const TString dummyName[1] = {""};
const int N_SigmaBarlow = 2;
#endif
