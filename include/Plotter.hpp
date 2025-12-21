#pragma once

#include <matplot/matplot.h>
#include <string>
#include <vector>
#include <algorithm>

template<typename T>
void universalPlot(const std::vector<std::vector<double>> &axisesY, const std::vector<T> &axisX, std::vector<std::string> plotNames, std::string labelX, std::string labelY, std::string title) {

    auto fig1 = matplot::figure();
    fig1->size(1200, 780);

    auto ax = matplot::gca();
    if(axisesY.size() > 1)
        matplot::hold(ax, true);
    
    int counter = 0;
    for(const auto axisY: axisesY) {
        matplot::plot(ax, axisX, axisY)->line_width(2).display_name(plotNames[counter]);
        counter++;
    }

    ax->xlabel(labelX);
    ax->ylabel(labelY);
    ax->title(title);
    ax->x_axis().label_font_size(18);
    ax->y_axis().label_font_size(18);
    ax->font_size(12);
    ax->title_font_size_multiplier(1.0);

    auto l = matplot::legend("show");
    l->font_size(16);
    l->location(matplot::legend::general_alignment::topright); 
    
    
    matplot::grid(true);
    matplot::save(title + ".png");
    matplot::show();
}


template<typename T>
void universalPlotVec(const std::vector<std::vector<double>> &axisesY, const std::vector<std::vector<T>> &axisesX, std::vector<std::string> plotNames, std::string labelX, std::string labelY, std::string title) {

    auto fig1 = matplot::figure();
    fig1->size(1200, 780);

    auto ax = matplot::gca();
    if(axisesY.size() > 1) matplot::hold(ax, true);
    
    if(axisesY.size() != axisesX.size()) return;
    
    for(int t = 0; t < axisesY.size(); t++) {
        matplot::plot(ax, axisesX[t], axisesY[t])->line_width(2).display_name(plotNames[t]);
    }

    ax->xlabel(labelX);
    ax->ylabel(labelY);
    ax->title(title);
    ax->x_axis().label_font_size(18);
    ax->y_axis().label_font_size(18);
    ax->font_size(12);
    ax->title_font_size_multiplier(1.0);

    auto l = matplot::legend("show");
    l->font_size(16);
    l->location(matplot::legend::general_alignment::topright); 
    
    
    matplot::grid(true);
    matplot::save(title + ".png");
    matplot::show();
}


// Plot with optional sigma (std-dev) shading for each series.
// sigmas: same shape as axisesY; if sigmas[i] is empty or has wrong size, no shading for that series.
template<typename T>
void universalPlotWithSigma(const std::vector<std::vector<double>> &axisesY, const std::vector<T> &axisX, const std::vector<std::vector<double>> &sigmas, std::vector<std::string> plotNames, std::string labelX, std::string labelY, std::string title) {

    auto fig1 = matplot::figure();
    fig1->size(1200, 780);

    auto ax = matplot::gca();
    if(axisesY.size() > 1)
        matplot::hold(ax, true);

    int counter = 0;
    for (size_t i = 0; i < axisesY.size(); ++i) {
        const auto &y = axisesY[i];
        // convert axisX elements to double
        std::vector<double> x;
        x.reserve(axisX.size());
        for (const auto &xx: axisX) x.push_back(static_cast<double>(xx));

        // draw shaded sigma band if available and sizes match
        if (i < sigmas.size() && sigmas[i].size() == y.size()) {
            std::vector<double> lower(y.size()), upper(y.size());
            for (size_t k = 0; k < y.size(); ++k) {
                lower[k] = y[k] - sigmas[i][k];
                upper[k] = y[k] + sigmas[i][k];
            }
            // build concatenated polygon coordinates: x + reverse(x), upper + reverse(lower)
            std::vector<double> xpoly;
            xpoly.reserve(x.size() * 2);
            xpoly.insert(xpoly.end(), x.begin(), x.end());
            for (auto it = x.rbegin(); it != x.rend(); ++it) xpoly.push_back(*it);

            std::vector<double> ypoly;
            ypoly.reserve(y.size() * 2);
            ypoly.insert(ypoly.end(), upper.begin(), upper.end());
            for (auto it = lower.rbegin(); it != lower.rend(); ++it) ypoly.push_back(*it);

            (void)matplot::fill(ax, xpoly, ypoly);
        }

        matplot::plot(ax, x, y)->line_width(2).display_name(plotNames[counter]);
        ++counter;
    }

    ax->xlabel(labelX);
    ax->ylabel(labelY);
    ax->title(title);
    ax->x_axis().label_font_size(18);
    ax->y_axis().label_font_size(18);
    ax->font_size(12);
    ax->title_font_size_multiplier(1.0);

    auto l = matplot::legend("show");
    l->font_size(16);
    l->location(matplot::legend::general_alignment::topright);

    matplot::grid(true);
    matplot::save(title + ".png");
    matplot::show();
}


template<typename T>
void universalPlotVecWithSigma(const std::vector<std::vector<double>> &axisesY, const std::vector<std::vector<T>> &axisesX, const std::vector<std::vector<double>> &sigmas, std::vector<std::string> plotNames, std::string labelX, std::string labelY, std::string title) {
    if (axisesY.size() != axisesX.size()) return;

    auto fig1 = matplot::figure();
    fig1->size(1200, 780);

    auto ax = matplot::gca();
    if(axisesY.size() > 1) matplot::hold(ax, true);

    for (size_t t = 0; t < axisesY.size(); t++) {
        const auto &y = axisesY[t];
        // convert x
        std::vector<double> x;
        x.reserve(axisesX[t].size());
        for (const auto &xx: axisesX[t]) x.push_back(static_cast<double>(xx));

        if (t < sigmas.size() && sigmas[t].size() == y.size()) {
            std::vector<double> lower(y.size()), upper(y.size());
            for (size_t k = 0; k < y.size(); ++k) {
                lower[k] = y[k] - sigmas[t][k];
                upper[k] = y[k] + sigmas[t][k];
            }
            std::vector<double> xpoly;
            xpoly.reserve(x.size() * 2);
            xpoly.insert(xpoly.end(), x.begin(), x.end());
            for (auto it = x.rbegin(); it != x.rend(); ++it) xpoly.push_back(*it);

            std::vector<double> ypoly;
            ypoly.reserve(y.size() * 2);
            ypoly.insert(ypoly.end(), upper.begin(), upper.end());
            for (auto it = lower.rbegin(); it != lower.rend(); ++it) ypoly.push_back(*it);

            (void)matplot::fill(ax, xpoly, ypoly);
        }

        matplot::plot(ax, x, y)->line_width(2).display_name(plotNames[t]);
    }

    ax->xlabel(labelX);
    ax->ylabel(labelY);
    ax->title(title);
    ax->x_axis().label_font_size(18);
    ax->y_axis().label_font_size(18);
    ax->font_size(12);
    ax->title_font_size_multiplier(1.0);

    auto l = matplot::legend("show");
    l->font_size(16);
    l->location(matplot::legend::general_alignment::topright); 
    
    matplot::grid(true);
    matplot::save(title + ".png");
    matplot::show();
}
