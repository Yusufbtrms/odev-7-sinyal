// Ödev 7 (görev 5 ve 6): frekans filtresi (sadece başlık dosyası).
// sinyal.hpp ile aynı klasörde olmalıdır.

#pragma once

#include "sinyal.hpp"

namespace sinyal {

// Frekans filtresi:
//   1) Sinyalin DFT'si alınır.
//   2) Genliği (harmonik vektöründeki değeri) eşiğin ALTINDA kalan harmonikler atılır (sıfırlanır).
//      Gerçek bir sinyalin DFT'si eşlenik simetriktir (X[N-k] = X[k]'nin eşleniği), bu yüzden k. harmoniği
//      atarken N-k karşılığı da atılır, korurken ikisi birlikte korunur; böylece sonuç gerçek sayılardan oluşur.
//   3) Kalan harmoniklerle ters DFT yapılır ve süzülmüş sinyal elde edilir.
//   threshold  : harmonik genliği eşiği (sinyalin genlik birimiyle aynı birim)
//   keptBins   : verilirse korunan harmonik numaraları (k = 0..N/2) döndürülür
//   keptAmp    : verilirse korunan harmoniklerin genlikleri döndürülür
inline std::vector<double> frequencyFilter(const std::vector<double>& x, double threshold, std::vector<int>* keptBins = nullptr,
                                           std::vector<double>* keptAmp = nullptr) {
    if (!isPowerOfTwo(x.size())) throw std::invalid_argument("vektor boyutu 2'nin kuvveti olmali");
    const size_t N = x.size(), half = N / 2;
    const std::vector<std::complex<double>> X = dftComplex(x);
    std::vector<std::complex<double>> Y(N, std::complex<double>(0.0, 0.0));
    if (keptBins) keptBins->clear();
    if (keptAmp) keptAmp->clear();
    for (size_t k = 0; k <= half; ++k) {
        const double mag = std::abs(X[k]) / static_cast<double>(N);
        const double amp = (k == 0 || k == half) ? mag : 2.0 * mag;
        if (amp >= threshold) {
            Y[k] = X[k];
            if (k != 0 && k != half) Y[N - k] = X[N - k];
            if (keptBins) keptBins->push_back(static_cast<int>(k));
            if (keptAmp) keptAmp->push_back(amp);
        }
    }
    return idftReal(Y);
}

inline std::vector<double> addSignals(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("vektor boyutlari ayni olmali");
    std::vector<double> s(a.size());
    for (size_t i = 0; i < a.size(); ++i) s[i] = a[i] + b[i];
    return s;
}

}  // namespace sinyal
