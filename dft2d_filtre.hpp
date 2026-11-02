// Ödev 7 (görev 9): frekans uzayında resim süzme (sadece başlık dosyası).
// dft2d.hpp ile aynı klasörde olmalıdır.
//
// Frekans filtresi: resmin 2 boyutlu DFT'si alınır, bazı harmonikler atılır (ya da küçültülür), ters DFT ile resim
// yeniden üretilir. Üç süzme yöntemi:
//   thresholdFilter2D : harmonik genliği eşiğin altında kalan harmonikler atılır (7. ödevin 6. görevindeki filtrenin resim hali)
//   lowPassFilter2D   : sıfır frekanstan uzaklığı (yarıçap) R'den büyük harmonikler atılır (ideal alçak geçiren)
//   gaussianLowPass2D : her harmonik, yarıçapına göre exp(-r^2 / (2 sigma^2)) ile çarpılır (yumuşak alçak geçiren)

#pragma once

#include "dft2d.hpp"

namespace sinyal {

inline double freqRadius(int u, int v, int W, int H) {
    const int fu = (u <= W / 2) ? u : u - W, fv = (v <= H / 2) ? v : v - H;
    return std::sqrt(static_cast<double>(fu) * fu + static_cast<double>(fv) * fv);
}

inline Spectrum2D thresholdFilter2D(const Spectrum2D& S, double threshold, long long* kept = nullptr) {
    Spectrum2D Y = S;
    long long cnt = 0;
    for (int v = 0; v < S.H; ++v)
        for (int u = 0; u < S.W; ++u) {
            const size_t i = static_cast<size_t>(v) * S.W + u;
            const double mag = std::abs(S.X[i]) / (static_cast<double>(S.H) * S.W);
            const bool selfConj = (u == 0 || u == S.W / 2) && (v == 0 || v == S.H / 2);
            const double amp = selfConj ? mag : 2.0 * mag;
            if ((u == 0 && v == 0) || amp >= threshold) {
                ++cnt;
            } else {
                Y.X[i] = std::complex<double>(0.0, 0.0);
            }
        }
    if (kept) *kept = cnt;
    return Y;
}

inline Spectrum2D lowPassFilter2D(const Spectrum2D& S, double radius, long long* kept = nullptr) {
    Spectrum2D Y = S;
    long long cnt = 0;
    for (int v = 0; v < S.H; ++v)
        for (int u = 0; u < S.W; ++u) {
            const size_t i = static_cast<size_t>(v) * S.W + u;
            if (freqRadius(u, v, S.W, S.H) <= radius) ++cnt;
            else Y.X[i] = std::complex<double>(0.0, 0.0);
        }
    if (kept) *kept = cnt;
    return Y;
}

inline Spectrum2D gaussianLowPass2D(const Spectrum2D& S, double sigmaF) {
    Spectrum2D Y = S;
    for (int v = 0; v < S.H; ++v)
        for (int u = 0; u < S.W; ++u) {
            const double r = freqRadius(u, v, S.W, S.H);
            Y.X[static_cast<size_t>(v) * S.W + u] *= std::exp(-r * r / (2.0 * sigmaF * sigmaF));
        }
    return Y;
}

// Ters DFT sonucunu (CV_64F) 0..255 aralığına sınırlayıp en yakın tam sayıya yuvarlar
inline cv::Mat toImage8(const cv::Mat& f64) {
    cv::Mat out(f64.size(), CV_8UC1);
    for (int y = 0; y < f64.rows; ++y)
        for (int x = 0; x < f64.cols; ++x) out.at<uchar>(y, x) = cv::saturate_cast<uchar>(f64.at<double>(y, x));
    return out;
}

}  // namespace sinyal
