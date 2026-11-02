// Ödev 7 - Görev 7: 2 boyutlu DFT. Fonksiyon bir resmi alır ve harmonik matrisini döndürür.
//
// Test, kendi harmoniklerini bildiğimiz sentetik desenlerle ve gerçek bir resimle yapılır:
//   A) dikey çizgiler (u = 10, v = 0):     iki tepe beklenir: (+10, 0) ve eşleniği (-10, 0), genlik 100
//   B) çapraz çizgiler (u = 8, v = 16):    tepe (+8, +16) ve eşleniği, genlik 100
//   C) iki frekansın toplamı:              (10, 0) genlik 60 ve (0, 24) genlik 30
//   D) beyaz dikdörtgen (64 x 32):         keskin kenarlı bir şeklin spektrumu (sinc benzeri yapı)
//   E) gerçek resim (256 x 256'ya küçültülmüş): doğal bir resmin spektrumu
// Her test için: ters DFT ile geri dönüşüm hatası, Parseval enerji kontrolü, en güçlü harmonikler,
// DC genliği (ortalama parlaklık), düşük frekanslardaki enerji oranı ve çalışma süresi yazılır.
// Her desen için [resim | log spektrum (sıfır frekans ortada)] görüntüsü kaydedilir.
//
// Kullanım: ./dft2d [cikti_dizini] [gercek_resim_yolu=images/Cat.jpg]
// Derleme:  g++ -O2 odev7_7_dft2d.cpp -o dft2d.exe -std=c++17 $(pkg-config --cflags --libs opencv5)
// (sinyal.hpp ve dft2d.hpp bu dosyayla aynı klasörde olmalıdır)

#include "dft2d.hpp"

#include <opencv2/imgcodecs.hpp>

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fsys = std::filesystem;

struct Component {
    int u0, v0;
    double amp;
};

static cv::Mat makePattern(int N, const std::vector<Component>& comps) {
    cv::Mat img(N, N, CV_8UC1);
    for (int y = 0; y < N; ++y)
        for (int x = 0; x < N; ++x) {
            double v = 128.0;
            for (const Component& c : comps) v += c.amp * std::cos(2.0 * sinyal::PI * (c.u0 * x + c.v0 * y) / N);
            img.at<uchar>(y, x) = cv::saturate_cast<uchar>(v);
        }
    return img;
}

static cv::Mat makeRect(int N, int w, int h) {
    cv::Mat img = cv::Mat::zeros(N, N, CV_8UC1);
    for (int y = (N - h) / 2; y < (N + h) / 2; ++y)
        for (int x = (N - w) / 2; x < (N + w) / 2; ++x) img.at<uchar>(y, x) = 255;
    return img;
}

static void analyse(const std::string& name, const std::string& file, const cv::Mat& img, const std::vector<Component>* expected, const fsys::path& outDir) {
    std::cout << "---- " << name << " (" << img.cols << "x" << img.rows << ") ----\n";
    const auto t0 = std::chrono::steady_clock::now();
    const sinyal::Spectrum2D S = sinyal::dft2dComplex(img);
    const auto t1 = std::chrono::steady_clock::now();
    const cv::Mat amp = sinyal::harmonicsFromSpectrum(S);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Sure (2 boyutlu DFT): " << std::setprecision(1) << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms\n" << std::setprecision(4);

    // Dogrulamalar
    const cv::Mat back = sinyal::idft2dReal(S);
    double maxErr = 0.0;
    cv::Mat f64;
    img.convertTo(f64, CV_64F);
    for (int y = 0; y < img.rows; ++y)
        for (int x = 0; x < img.cols; ++x) maxErr = std::max(maxErr, std::fabs(back.at<double>(y, x) - f64.at<double>(y, x)));
    double eSpace = 0.0, eFreq = 0.0;
    for (int y = 0; y < img.rows; ++y)
        for (int x = 0; x < img.cols; ++x) eSpace += f64.at<double>(y, x) * f64.at<double>(y, x);
    for (const auto& z : S.X) eFreq += std::norm(z);
    eFreq /= static_cast<double>(img.rows) * img.cols;
    std::cout << "  Ters DFT ile geri donusum: en buyuk fark = " << std::scientific << maxErr << std::fixed << "\n";
    std::cout << "  Parseval: sum|f|^2 = " << std::setprecision(1) << eSpace << ", sum|X|^2/(H*W) = " << eFreq << std::setprecision(4) << " (esit olmali)\n";
    std::cout << "  DC genligi = " << amp.at<double>(0, 0) << " (resmin ortalama parlakligi: " << cv::mean(img)[0] << ")\n";

    // Dusuk frekans enerji orani (|fu| <= 16 ve |fv| <= 16, DC haric)
    double eLow = 0.0, eAll = 0.0;
    for (int v = 0; v < S.H; ++v)
        for (int u = 0; u < S.W; ++u) {
            if (u == 0 && v == 0) continue;
            const double e = std::norm(S.X[static_cast<size_t>(v) * S.W + u]);
            eAll += e;
            const int fu = (u <= S.W / 2) ? u : u - S.W, fv = (v <= S.H / 2) ? v : v - S.H;
            if (std::abs(fu) <= 16 && std::abs(fv) <= 16) eLow += e;
        }
    std::cout << "  Enerjinin (DC haric) |fu|,|fv| <= 16 dusuk frekans bolgesindeki orani = " << std::setprecision(2) << 100.0 * eLow / eAll << " %\n" << std::setprecision(4);

    std::cout << "  En guclu 4 harmonik (DC haric, eslenik ciftlerden biri; genlik = sinus genligi):";
    for (const sinyal::SpecPeak& p : sinyal::topPeaks(amp, 4)) std::cout << " (fu=" << p.fu << ", fv=" << p.fv << ", " << std::setprecision(3) << p.amp << ")" << std::setprecision(4);
    std::cout << "\n";
    if (expected) {
        std::cout << "  Beklenen:";
        for (const Component& c : *expected) std::cout << " (fu=" << c.u0 << ", fv=" << c.v0 << ", " << c.amp << ")";
        std::cout << "\n";
    }

    // Gorsel: [resim | log spektrum], her biri 2 kat buyutulmus
    cv::Mat left, right, both;
    cv::resize(img, left, cv::Size(), 2.0, 2.0, cv::INTER_NEAREST);
    cv::resize(sinyal::spectrumToImage(amp), right, cv::Size(), 2.0, 2.0, cv::INTER_NEAREST);
    cv::Mat l3, r3;
    cv::cvtColor(left, l3, cv::COLOR_GRAY2BGR);
    cv::cvtColor(right, r3, cv::COLOR_GRAY2BGR);
    cv::putText(l3, "goruntu", cv::Point(8, 22), cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(0, 255, 255), 1, cv::LINE_AA);
    cv::putText(r3, "log spektrum (DC ortada)", cv::Point(8, 22), cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(0, 255, 255), 1, cv::LINE_AA);
    cv::hconcat(std::vector<cv::Mat>{l3, r3}, both);
    cv::imwrite((outDir / file).string(), both);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    const fsys::path outDir = (argc >= 2) ? fsys::path(argv[1]) : fsys::path("cikti_dft2d");
    const std::string realPath = (argc >= 3) ? argv[2] : "images/Cat.jpg";
    fsys::create_directories(outDir);
    const int N = 256;

    std::cout << "2 BOYUTLU DFT (harmonik matrisi; hucre (fu, fv) = o frekanstaki sinusun genligi, fu: yatay, fv: dikey frekans)\n\n";

    const std::vector<Component> a = {{10, 0, 100.0}};
    analyse("A) dikey cizgiler", "7_a_dikey_cizgiler.png", makePattern(N, a), &a, outDir);
    const std::vector<Component> b = {{8, 16, 100.0}};
    analyse("B) capraz cizgiler", "7_b_capraz_cizgiler.png", makePattern(N, b), &b, outDir);
    const std::vector<Component> c = {{10, 0, 60.0}, {0, 24, 30.0}};
    analyse("C) iki frekansin toplami", "7_c_iki_frekans.png", makePattern(N, c), &c, outDir);
    analyse("D) beyaz dikdortgen 64x32", "7_d_dikdortgen.png", makeRect(N, 64, 32), nullptr, outDir);

    cv::Mat real = cv::imread(realPath, cv::IMREAD_GRAYSCALE);
    if (real.empty()) {
        std::cout << "(gercek resim okunamadi: " << realPath << ", E testi atlandi)\n";
    } else {
        cv::Mat small;
        cv::resize(real, small, cv::Size(N, N), 0, 0, cv::INTER_AREA);
        analyse("E) gercek resim (" + realPath + ", 256x256'ya kucultuldu)", "7_e_gercek_resim.png", small, nullptr, outDir);
    }

    std::cout << "Kaydedilenler (" << outDir.string() << "):\n"
              << "  7_a_dikey_cizgiler.png, 7_b_capraz_cizgiler.png, 7_c_iki_frekans.png, 7_d_dikdortgen.png, 7_e_gercek_resim.png\n"
              << "  (her biri: solda resim, sagda log spektrum; spektrumda sifir frekans ortada)\n";
    return 0;
}
