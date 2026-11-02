// Ödev 7 - Sinyal işleme fonksiyonları (sadece başlık dosyası).
// Bu dosyayı programla AYNI klasöre koyun; ayrıca derlemenize gerek yoktur.
//
// İçindekiler:
//   Madde 1: fillSine()          : faz, genlik, frekans (ya da periyot) ile bir vektörü sinüs değerleriyle doldurur
//   Madde 2: addNoise()          : genliğe ve gürültü oranına göre vektöre rastgele gürültü ekler
//   Madde 3: gaussianFilter(), medianFilter() : tek boyutlu Gauss ve medyan filtre
//   Madde 4: dftHarmonics()      : DFT ile frekans uzayına geçip harmonik (genlik) vektörünü döndürür
//   Yardımcılar: addImpulseNoise(), idftReal(), plotSeries() (OpenCV ile grafik çizimi)

#pragma once

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sinyal {

const double PI = 3.14159265358979323846;

inline bool isPowerOfTwo(size_t n) { return n > 0 && (n & (n - 1)) == 0; }

// ---------------------------------------------------------------------------------------------
// Madde 1: sinüs üretici
// v[n] = genlik * sin(2*pi*frekans*n/ornekleme + faz)
//   v           : doldurulacak vektör, boyutu 2'nin kuvveti olmalıdır
//   amplitude   : genlik
//   phaseRad    : faz (radyan)
//   frequencyHz : frekans (Hz)
//   sampleRateHz: örnekleme frekansı (Hz); vektörün her elemanı 1/sampleRateHz saniye aralıklı bir örnektir
// ---------------------------------------------------------------------------------------------
inline void fillSine(std::vector<double>& v, double amplitude, double phaseRad, double frequencyHz, double sampleRateHz) {
    if (!isPowerOfTwo(v.size())) throw std::invalid_argument("vektor boyutu 2'nin kuvveti olmali");
    for (size_t n = 0; n < v.size(); ++n) v[n] = amplitude * std::sin(2.0 * PI * frequencyHz * static_cast<double>(n) / sampleRateHz + phaseRad);
}

// Aynı işlev, periyot (örnek sayısı cinsinden) ile: periyot = 1 / frekans (frekans "örnek başına döngü")
inline void fillSineByPeriod(std::vector<double>& v, double amplitude, double phaseRad, double periodInSamples) {
    fillSine(v, amplitude, phaseRad, 1.0 / periodInSamples, 1.0);
}

// ---------------------------------------------------------------------------------------------
// Madde 2: gürültü ekleme
// Vektörün genliği (en büyük mutlak değer) A ise, her elemana  oran * A * u  eklenir;
// u, [-1, 1] aralığında düzgün dağılımlı rastgele bir sayıdır. Yani en büyük gürültü oran * A kadardır.
// ---------------------------------------------------------------------------------------------
inline std::vector<double> addNoise(const std::vector<double>& v, double ratio, unsigned seed = 42) {
    if (ratio < 0.0 || ratio > 1.0) throw std::invalid_argument("gurultu orani [0, 1] araliginda olmali");
    double amplitude = 0.0;
    for (double x : v) amplitude = std::max(amplitude, std::fabs(x));
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<double> out(v.size());
    for (size_t i = 0; i < v.size(); ++i) out[i] = v[i] + ratio * amplitude * u(rng);
    return out;
}

// Ek deney: dürtü (impuls) gürültüsü. Elemanların "fraction" oranı rastgele seçilir ve +-1.5*genlik değerine ayarlanır.
inline std::vector<double> addImpulseNoise(const std::vector<double>& v, double fraction, unsigned seed = 7) {
    double amplitude = 0.0;
    for (double x : v) amplitude = std::max(amplitude, std::fabs(x));
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<double> out = v;
    for (size_t i = 0; i < v.size(); ++i)
        if (u(rng) < fraction) out[i] = (u(rng) < 0.5 ? -1.5 : 1.5) * amplitude;
    return out;
}

// ---------------------------------------------------------------------------------------------
// Madde 3: filtreler (kenarlarda yansıtma)
// ---------------------------------------------------------------------------------------------
inline int reflectIndex(int i, int n) {
    if (n == 1) return 0;
    while (i < 0 || i >= n) {
        if (i < 0) i = -i;
        if (i >= n) i = 2 * (n - 1) - i;
    }
    return i;
}

// Gauss filtresi: ağırlıklar exp(-d^2 / (2 sigma^2)), toplamı 1'e normalize; size tek sayı olmalıdır
inline std::vector<double> gaussianFilter(const std::vector<double>& v, int size, double sigma) {
    if (size % 2 == 0 || size < 1) throw std::invalid_argument("Gauss pencere boyutu tek sayi olmali");
    const int r = size / 2, n = static_cast<int>(v.size());
    std::vector<double> k(size);
    double sum = 0.0;
    for (int j = 0; j < size; ++j) {
        const double d = j - r;
        k[j] = std::exp(-d * d / (2.0 * sigma * sigma));
        sum += k[j];
    }
    for (double& w : k) w /= sum;
    std::vector<double> out(v.size(), 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < size; ++j) s += k[j] * v[reflectIndex(i + j - r, n)];
        out[i] = s;
    }
    return out;
}

// Medyan filtresi: penceredeki değerler sıralanır, ortadaki değer alınır; size tek sayı olmalıdır
inline std::vector<double> medianFilter(const std::vector<double>& v, int size) {
    if (size % 2 == 0 || size < 1) throw std::invalid_argument("medyan pencere boyutu tek sayi olmali");
    const int r = size / 2, n = static_cast<int>(v.size());
    std::vector<double> out(v.size(), 0.0), window(size);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < size; ++j) window[j] = v[reflectIndex(i + j - r, n)];
        std::sort(window.begin(), window.end());
        out[i] = window[r];
    }
    return out;
}

// ---------------------------------------------------------------------------------------------
// Madde 4: DFT
//   X[k] = sum_n x[n] * exp(-j*2*pi*k*n/N)      (doğrudan tanım, O(N^2); sin/cos tablosu bir kere hesaplanır)
// ---------------------------------------------------------------------------------------------
inline std::vector<std::complex<double>> dftComplex(const std::vector<double>& x) {
    const size_t N = x.size();
    std::vector<double> cosT(N), sinT(N);
    for (size_t m = 0; m < N; ++m) {
        cosT[m] = std::cos(2.0 * PI * static_cast<double>(m) / static_cast<double>(N));
        sinT[m] = std::sin(2.0 * PI * static_cast<double>(m) / static_cast<double>(N));
    }
    std::vector<std::complex<double>> X(N);
    for (size_t k = 0; k < N; ++k) {
        double re = 0.0, im = 0.0;
        for (size_t n = 0; n < N; ++n) {
            const size_t idx = (k * n) % N;
            re += x[n] * cosT[idx];
            im -= x[n] * sinT[idx];
        }
        X[k] = std::complex<double>(re, im);
    }
    return X;
}

// Ters DFT (doğrulama için): x[n] = (1/N) * sum_k X[k] * exp(+j*2*pi*k*n/N)
inline std::vector<double> idftReal(const std::vector<std::complex<double>>& X) {
    const size_t N = X.size();
    std::vector<double> cosT(N), sinT(N);
    for (size_t m = 0; m < N; ++m) {
        cosT[m] = std::cos(2.0 * PI * static_cast<double>(m) / static_cast<double>(N));
        sinT[m] = std::sin(2.0 * PI * static_cast<double>(m) / static_cast<double>(N));
    }
    std::vector<double> x(N, 0.0);
    for (size_t n = 0; n < N; ++n) {
        double s = 0.0;
        for (size_t k = 0; k < N; ++k) {
            const size_t idx = (k * n) % N;
            s += X[k].real() * cosT[idx] - X[k].imag() * sinT[idx];
        }
        x[n] = s / static_cast<double>(N);
    }
    return x;
}

// Harmonik vektörü: k = 0..N/2 için her harmoniğin GENLİĞİ (sinüs genliği cinsinden).
//   amp[0] = |X[0]| / N (DC bileşeni), amp[N/2] = |X[N/2]| / N, diğerleri 2*|X[k]| / N.
// phaseOut verilirse her harmoniğin DFT fazı (radyan) da döndürülür.
inline std::vector<double> dftHarmonics(const std::vector<double>& x, std::vector<double>* phaseOut = nullptr) {
    if (!isPowerOfTwo(x.size())) throw std::invalid_argument("vektor boyutu 2'nin kuvveti olmali");
    const size_t N = x.size(), half = N / 2;
    const std::vector<std::complex<double>> X = dftComplex(x);
    std::vector<double> amp(half + 1);
    if (phaseOut) phaseOut->assign(half + 1, 0.0);
    for (size_t k = 0; k <= half; ++k) {
        const double mag = std::abs(X[k]) / static_cast<double>(N);
        amp[k] = (k == 0 || k == half) ? mag : 2.0 * mag;
        if (phaseOut) (*phaseOut)[k] = std::arg(X[k]);
    }
    return amp;
}

// ---------------------------------------------------------------------------------------------
// Grafik çizimi (OpenCV ile; yazılar ASCII)
// ---------------------------------------------------------------------------------------------
struct Series {
    std::string name;
    std::vector<double> x;  // boşsa 0, 1, 2, ... indeksleri kullanılır
    std::vector<double> y;
    cv::Scalar color;
    bool line = true;
    bool dots = false;
    bool stems = false;  // 0 çizgisinden değere dikey çubuk
    int thickness = 1;
};

inline std::string fmt(double v, int prec) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(prec) << v;
    return o.str();
}

inline cv::Mat plotSeries(const std::vector<Series>& series, const std::string& title, double yMin, double yMax, double xMin, double xMax,
                          const std::string& xLabel, int width = 1100, int height = 420) {
    cv::Mat img(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
    const int left = 62, right = 20, top = 52, bottom = 46;
    const int pw = width - left - right, ph = height - top - bottom;
    auto PX = [&](double x) { return left + static_cast<int>(std::lround((x - xMin) / (xMax - xMin) * pw)); };
    auto PY = [&](double y) { return top + static_cast<int>(std::lround((yMax - y) / (yMax - yMin) * ph)); };
    const cv::Scalar gridC(225, 225, 225), textC(70, 70, 70);

    for (int i = 0; i <= 4; ++i) {
        const double yv = yMin + (yMax - yMin) * i / 4.0;
        const int py = PY(yv);
        cv::line(img, cv::Point(left, py), cv::Point(left + pw, py), gridC, 1);
        cv::putText(img, fmt(yv, 2), cv::Point(6, py + 4), cv::FONT_HERSHEY_SIMPLEX, 0.42, textC, 1, cv::LINE_AA);
    }
    for (int i = 0; i <= 8; ++i) {
        const double xv = xMin + (xMax - xMin) * i / 8.0;
        const int px = PX(xv);
        cv::line(img, cv::Point(px, top), cv::Point(px, top + ph), gridC, 1);
        cv::putText(img, fmt(xv, 0), cv::Point(px - 10, top + ph + 16), cv::FONT_HERSHEY_SIMPLEX, 0.42, textC, 1, cv::LINE_AA);
    }
    if (yMin < 0.0 && yMax > 0.0) cv::line(img, cv::Point(left, PY(0.0)), cv::Point(left + pw, PY(0.0)), cv::Scalar(170, 170, 170), 1);
    cv::rectangle(img, cv::Rect(left, top, pw, ph), cv::Scalar(120, 120, 120), 1);
    cv::putText(img, title, cv::Point(left, 20), cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(30, 30, 30), 1, cv::LINE_AA);
    cv::putText(img, xLabel, cv::Point(left + pw / 2 - 50, height - 10), cv::FONT_HERSHEY_SIMPLEX, 0.45, textC, 1, cv::LINE_AA);

    for (const Series& s : series) {
        cv::Point prev;
        bool hasPrev = false;
        for (size_t i = 0; i < s.y.size(); ++i) {
            const double xv = s.x.empty() ? static_cast<double>(i) : s.x[i];
            if (xv < xMin || xv > xMax) {
                hasPrev = false;
                continue;
            }
            const cv::Point p(PX(xv), PY(s.y[i]));
            if (s.stems) cv::line(img, cv::Point(p.x, PY(0.0)), p, s.color, 2, cv::LINE_AA);
            if (s.line && hasPrev) cv::line(img, prev, p, s.color, s.thickness, cv::LINE_AA);
            if (s.dots) cv::circle(img, p, 3, s.color, -1, cv::LINE_AA);
            prev = p;
            hasPrev = true;
        }
    }

    int lx = left;  // aciklama (legend)
    for (const Series& s : series) {
        int baseline = 0;
        const cv::Size ts = cv::getTextSize(s.name, cv::FONT_HERSHEY_SIMPLEX, 0.45, 1, &baseline);
        cv::line(img, cv::Point(lx, 38), cv::Point(lx + 22, 38), s.color, 3, cv::LINE_AA);
        cv::putText(img, s.name, cv::Point(lx + 28, 43), cv::FONT_HERSHEY_SIMPLEX, 0.45, textC, 1, cv::LINE_AA);
        lx += 28 + ts.width + 24;
    }
    return img;
}

}  // namespace sinyal
