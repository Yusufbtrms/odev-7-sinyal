// Ödev 7 (görev 7): 2 boyutlu DFT (sadece başlık dosyası).
// sinyal.hpp ile aynı klasörde olmalıdır.
//
// 2 boyutlu DFT, tanımı gereği AYRILABİLİR: önce her satırın tek boyutlu DFT'si alınır, sonra çıkan sonucun
// her sütununun tek boyutlu DFT'si alınır.
//   X[v][u] = sum_y sum_x f[y][x] * exp(-j*2*pi*(u*x/W + v*y/H))
//   u: yatay frekans (sütun), v: dikey frekans (satır)
// Resmin en ve boyu 2'nin kuvveti olmalıdır (vektörlerdeki koşulla aynı).

#pragma once

#include "sinyal.hpp"

#include <opencv2/core.hpp>

namespace sinyal {

struct Spectrum2D {
    int H = 0, W = 0;
    std::vector<std::complex<double>> X;  // X[v * W + u]
};

// Karmaşık 2 boyutlu DFT (doğrudan tanım, her eksen için O(N^2); sin/cos tabloları bir kere hesaplanır)
inline Spectrum2D dft2dComplex(const cv::Mat& img) {
    if (img.channels() != 1) throw std::invalid_argument("tek kanalli (gri) resim gerekli");
    const int H = img.rows, W = img.cols;
    if (!isPowerOfTwo(static_cast<size_t>(H)) || !isPowerOfTwo(static_cast<size_t>(W))) throw std::invalid_argument("resmin eni ve boyu 2'nin kuvveti olmali");

    cv::Mat f;
    img.convertTo(f, CV_64F);

    std::vector<double> cosW(W), sinW(W), cosH(H), sinH(H);
    for (int m = 0; m < W; ++m) {
        cosW[m] = std::cos(2.0 * PI * m / W);
        sinW[m] = std::sin(2.0 * PI * m / W);
    }
    for (int m = 0; m < H; ++m) {
        cosH[m] = std::cos(2.0 * PI * m / H);
        sinH[m] = std::sin(2.0 * PI * m / H);
    }

    // Adım 1: her satırın DFT'si (girdi gerçek sayı)
    std::vector<double> re1(static_cast<size_t>(H) * W), im1(static_cast<size_t>(H) * W);
    for (int y = 0; y < H; ++y) {
        const double* row = f.ptr<double>(y);
        for (int u = 0; u < W; ++u) {
            double sr = 0.0, si = 0.0;
            for (int x = 0; x < W; ++x) {
                const int idx = (u * x) % W;
                sr += row[x] * cosW[idx];
                si -= row[x] * sinW[idx];
            }
            re1[static_cast<size_t>(y) * W + u] = sr;
            im1[static_cast<size_t>(y) * W + u] = si;
        }
    }

    // Adım 2: Adım 1 sonucunun her sütununun DFT'si (girdi karmaşık sayı)
    Spectrum2D S;
    S.H = H;
    S.W = W;
    S.X.assign(static_cast<size_t>(H) * W, std::complex<double>(0.0, 0.0));
    std::vector<double> colRe(H), colIm(H);
    for (int u = 0; u < W; ++u) {
        for (int y = 0; y < H; ++y) {
            colRe[y] = re1[static_cast<size_t>(y) * W + u];
            colIm[y] = im1[static_cast<size_t>(y) * W + u];
        }
        for (int v = 0; v < H; ++v) {
            double sr = 0.0, si = 0.0;
            for (int y = 0; y < H; ++y) {
                const int idx = (v * y) % H;
                const double c = cosH[idx], s = sinH[idx];
                sr += colRe[y] * c + colIm[y] * s;  // (a + jb)(c - js) gerçel kısmı
                si += colIm[y] * c - colRe[y] * s;  // (a + jb)(c - js) sanal kısmı
            }
            S.X[static_cast<size_t>(v) * W + u] = std::complex<double>(sr, si);
        }
    }
    return S;
}

// Ters 2 boyutlu DFT (doğrulama için): f[y][x] = (1/(H*W)) * sum_v sum_u X[v][u] * exp(+j*2*pi*(u*x/W + v*y/H))
inline cv::Mat idft2dReal(const Spectrum2D& S) {
    const int H = S.H, W = S.W;
    std::vector<double> cosW(W), sinW(W), cosH(H), sinH(H);
    for (int m = 0; m < W; ++m) {
        cosW[m] = std::cos(2.0 * PI * m / W);
        sinW[m] = std::sin(2.0 * PI * m / W);
    }
    for (int m = 0; m < H; ++m) {
        cosH[m] = std::cos(2.0 * PI * m / H);
        sinH[m] = std::sin(2.0 * PI * m / H);
    }
    // sütunlar boyunca ters DFT
    std::vector<double> re2(static_cast<size_t>(H) * W), im2(static_cast<size_t>(H) * W);
    for (int u = 0; u < W; ++u) {
        for (int y = 0; y < H; ++y) {
            double sr = 0.0, si = 0.0;
            for (int v = 0; v < H; ++v) {
                const int idx = (v * y) % H;
                const std::complex<double>& z = S.X[static_cast<size_t>(v) * W + u];
                sr += z.real() * cosH[idx] - z.imag() * sinH[idx];
                si += z.real() * sinH[idx] + z.imag() * cosH[idx];
            }
            re2[static_cast<size_t>(y) * W + u] = sr;
            im2[static_cast<size_t>(y) * W + u] = si;
        }
    }
    // satırlar boyunca ters DFT (sadece gerçel kısım gerekir)
    cv::Mat out(H, W, CV_64F);
    for (int y = 0; y < H; ++y) {
        double* row = out.ptr<double>(y);
        for (int x = 0; x < W; ++x) {
            double s = 0.0;
            for (int u = 0; u < W; ++u) {
                const int idx = (u * x) % W;
                s += re2[static_cast<size_t>(y) * W + u] * cosW[idx] - im2[static_cast<size_t>(y) * W + u] * sinW[idx];
            }
            row[x] = s / (static_cast<double>(H) * W);
        }
    }
    return out;
}

// Spektrumdan harmonik matrisi: H x W, CV_64F.
//   Her hücre, o frekanstaki sinüsün GENLİĞİNİ verir (tek boyutlu dftHarmonics ile aynı anlam):
//     DC (0,0) ve kendi eşleniği olan hücreler (u ve v, 0 ya da yarım boyut): |X| / (H*W)
//     diğer hücreler: 2 * |X| / (H*W); eşlenik hücre (-u, -v) aynı değeri taşır.
//   DC hücresi resmin ortalama parlaklığıdır.
// phaseOut verilirse her hücrenin DFT fazı (radyan) da döndürülür.
inline cv::Mat harmonicsFromSpectrum(const Spectrum2D& S, cv::Mat* phaseOut = nullptr) {
    const int H = S.H, W = S.W;
    cv::Mat amp(H, W, CV_64F);
    if (phaseOut) phaseOut->create(H, W, CV_64F);
    for (int v = 0; v < H; ++v) {
        for (int u = 0; u < W; ++u) {
            const std::complex<double>& z = S.X[static_cast<size_t>(v) * W + u];
            const double mag = std::abs(z) / (static_cast<double>(H) * W);
            const bool selfConj = (u == 0 || u == W / 2) && (v == 0 || v == H / 2);
            amp.at<double>(v, u) = selfConj ? mag : 2.0 * mag;
            if (phaseOut) phaseOut->at<double>(v, u) = std::arg(z);
        }
    }
    return amp;
}

// 2 boyutlu DFT fonksiyonu: bir resmi alır, harmonik matrisini (frekans uzayı genlikleri) döndürür.
inline cv::Mat dft2dHarmonics(const cv::Mat& img, cv::Mat* phaseOut = nullptr) { return harmonicsFromSpectrum(dft2dComplex(img), phaseOut); }

// Gösterim: sıfır frekans (DC) merkeze kaydırılır, genlikler logaritmik ölçekle 0..255'e getirilir.
// DC hücresi diğerlerinden çok büyük olduğu için ölçek, DC dışındaki en büyük genliğe göre belirlenir.
inline cv::Mat spectrumToImage(const cv::Mat& amp) {
    const int H = amp.rows, W = amp.cols;
    double maxRef = 0.0;
    for (int v = 0; v < H; ++v)
        for (int u = 0; u < W; ++u)
            if (!(u == 0 && v == 0)) maxRef = std::max(maxRef, amp.at<double>(v, u));
    cv::Mat img(H, W, CV_8UC1);
    const double denom = std::log(1.0 + maxRef);
    for (int v = 0; v < H; ++v) {
        for (int u = 0; u < W; ++u) {
            const double val = denom > 0 ? std::log(1.0 + amp.at<double>(v, u)) / denom : 0.0;
            img.at<uchar>((v + H / 2) % H, (u + W / 2) % W) = cv::saturate_cast<uchar>(val * 255.0);
        }
    }
    return img;
}

struct SpecPeak {
    int fu, fv;  // merkezli frekans koordinatları (negatif olabilir)
    double amp;
};

// DC dışındaki en güçlü count tepeyi döndürür; eşlenik çiftlerden yalnızca biri listelenir
inline std::vector<SpecPeak> topPeaks(const cv::Mat& amp, int count) {
    const int H = amp.rows, W = amp.cols;
    std::vector<SpecPeak> all;
    for (int v = 0; v < H; ++v) {
        for (int u = 0; u < W; ++u) {
            if (u == 0 && v == 0) continue;
            const int fu = (u <= W / 2) ? u : u - W, fv = (v <= H / 2) ? v : v - H;
            const bool representative = (fv > 0) || (fv == 0 && fu > 0);
            const bool nyquistSelf = (u == W / 2 || u == 0) && (v == H / 2 || v == 0);  // kendi eşleniği
            if (representative || nyquistSelf) all.push_back({fu, fv, amp.at<double>(v, u)});
        }
    }
    std::sort(all.begin(), all.end(), [](const SpecPeak& a, const SpecPeak& b) { return a.amp > b.amp; });
    if (static_cast<int>(all.size()) > count) all.resize(count);
    return all;
}

}  // namespace sinyal
