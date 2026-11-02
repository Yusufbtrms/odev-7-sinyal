// Ödev 7: Sinüs üretme, gürültü ekleme, Gauss/medyan filtreleme ve DFT (madde 1-4).
//
//   Madde 1: N elemanlı (2'nin kuvveti) bir vektörü sinüs değerleriyle doldur, vektörü ve sinüs dalgasını çiz.
//   Madde 2: Vektöre genliğe ve gürültü oranına bağlı rastgele gürültü ekle, çiz.
//   Madde 3: Gürültülü veriyi Gauss ve medyan filtreyle temizle, çiz, ölç.
//   Madde 4: DFT ile harmonik vektörünü hesapla; temiz, gürültülü ve filtrelenmiş vektörler için karşılaştır.
//
// Varsayılan sinyal: N = 256 örnek, örnekleme frekansı 256 Hz (yani vektör tam 1 saniye),
// 8 Hz sinüs, genlik 2.0, faz 30 derece. Böylece 8 tam periyot vektöre sığar ve DFT'nin frekans
// çözünürlüğü 1 Hz olur (harmonik k, k Hz'e karşılık gelir).
//
// Kullanım: ./sinyal [cikti_dizini] [N=256] [gurultu_orani=0.3] [frekans_Hz=8] [genlik=2.0] [faz_derece=30]
// Derleme:  g++ -O2 odev7_sinyal.cpp -o sinyal.exe -std=c++17 $(pkg-config --cflags --libs opencv5)
// (sinyal.hpp bu dosyayla aynı klasörde olmalıdır)

#include "sinyal.hpp"

#include <opencv2/imgcodecs.hpp>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using sinyal::Series;

static double mse(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) s += (a[i] - b[i]) * (a[i] - b[i]);
    return s / static_cast<double>(a.size());
}

static double snrDb(const std::vector<double>& clean, const std::vector<double>& est) {
    double sig = 0.0, err = 0.0;
    for (size_t i = 0; i < clean.size(); ++i) {
        sig += clean[i] * clean[i];
        err += (est[i] - clean[i]) * (est[i] - clean[i]);
    }
    return err > 0.0 ? 10.0 * std::log10(sig / err) : 999.0;
}

struct SpecInfo {
    int peakBin = 0;
    double peakAmp = 0.0, phaseDeg = 0.0, floorMean = 0.0, floorMax = 0.0;
};

// Harmonik vektöründen tepe, faz tahmini ve gürültü tabanı (tepe dışındaki harmoniklerin ortalaması)
static SpecInfo analyse(const std::vector<double>& x, std::vector<double>* ampOut = nullptr) {
    std::vector<double> phase;
    const std::vector<double> amp = sinyal::dftHarmonics(x, &phase);
    SpecInfo s;
    for (size_t k = 1; k < amp.size(); ++k)
        if (amp[k] > s.peakAmp) {
            s.peakAmp = amp[k];
            s.peakBin = static_cast<int>(k);
        }
    // sin(w n + faz) icin DFT fazi = faz - 90 derece
    double ph = phase[s.peakBin] + sinyal::PI / 2.0;
    while (ph > sinyal::PI) ph -= 2.0 * sinyal::PI;
    while (ph < -sinyal::PI) ph += 2.0 * sinyal::PI;
    s.phaseDeg = ph * 180.0 / sinyal::PI;
    double sum = 0.0;
    int cnt = 0;
    for (size_t k = 1; k < amp.size(); ++k) {
        if (static_cast<int>(k) == s.peakBin) continue;
        sum += amp[k];
        s.floorMax = std::max(s.floorMax, amp[k]);
        ++cnt;
    }
    s.floorMean = cnt ? sum / cnt : 0.0;
    if (ampOut) *ampOut = amp;
    return s;
}

static double ampAtBin(const std::vector<double>& x, int bin) { return sinyal::dftHarmonics(x)[bin]; }

int main(int argc, char** argv) {
    const fs::path outDir = (argc >= 2) ? fs::path(argv[1]) : fs::path("cikti_sinyal");
    const int N = (argc >= 3) ? std::stoi(argv[2]) : 256;
    const double ratio = (argc >= 4) ? std::stod(argv[3]) : 0.3;
    const double freq = (argc >= 5) ? std::stod(argv[4]) : 8.0;
    const double amp = (argc >= 6) ? std::stod(argv[5]) : 2.0;
    const double phaseDeg = (argc >= 7) ? std::stod(argv[6]) : 30.0;
    fs::create_directories(outDir);

    if (!sinyal::isPowerOfTwo(static_cast<size_t>(N))) {
        std::cerr << "Hata: N = " << N << " 2'nin kuvveti degil (ornek: 64, 128, 256, 512)\n";
        return 1;
    }
    const double fs = static_cast<double>(N);  // 1 saniyelik pencere
    const double phase = phaseDeg * sinyal::PI / 180.0;
    const int bin = static_cast<int>(std::lround(freq * N / fs));
    const double yR = amp * 1.15 + 0.2;

    // =========================================================================================
    // MADDE 1: sinus vektoru
    // =========================================================================================
    std::vector<double> clean(N);
    sinyal::fillSine(clean, amp, phase, freq, fs);

    std::cout << "MADDE 1 - SINUS VEKTORU\n";
    std::cout << "  N = " << N << " (2'nin kuvveti), ornekleme = " << fs << " Hz, frekans = " << freq << " Hz (periyot = " << fs / freq << " ornek), genlik = " << amp
              << ", faz = " << phaseDeg << " derece\n";
    std::cout << "  Ilk 8 eleman:";
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < 8 && i < N; ++i) std::cout << " " << clean[i];
    std::cout << "\n  En buyuk = " << *std::max_element(clean.begin(), clean.end()) << ", en kucuk = " << *std::min_element(clean.begin(), clean.end()) << "\n\n";

    std::vector<double> xc, yc;  // surekli sinus dalgasi (cizim icin sik ornekli)
    for (double x = 0.0; x <= N - 1; x += 0.25) {
        xc.push_back(x);
        yc.push_back(amp * std::sin(2.0 * sinyal::PI * freq * x / fs + phase));
    }
    const cv::Scalar blue(200, 120, 0), red(0, 0, 220), green(0, 150, 0), gray(170, 170, 170), orange(0, 140, 255);
    {
        const Series curve{"sinus dalgasi (surekli)", xc, yc, blue, true, false, false, 2};
        const Series samples{"vektor elemanlari", {}, clean, red, false, true, false, 1};
        const Series samplesStem{"vektor elemanlari", {}, clean, red, false, true, true, 1};
        cv::Mat p1 = sinyal::plotSeries({curve, samples}, "Madde 1: sinus vektoru (tum vektor, N = " + std::to_string(N) + ")", -yR, yR, -1, N, "eleman no (n)");
        cv::Mat p2 = sinyal::plotSeries({curve, samplesStem}, "Madde 1: ilk 64 eleman ayrintili", -yR, yR, -1, 64, "eleman no (n)");
        cv::Mat both;
        cv::vconcat(std::vector<cv::Mat>{p1, p2}, both);
        cv::imwrite((outDir / "1_sinus_vektoru.png").string(), both);
    }

    // =========================================================================================
    // MADDE 2: gurultu ekleme
    // =========================================================================================
    const std::vector<double> noisy = sinyal::addNoise(clean, ratio, 42);
    double noiseMax = 0.0, noiseSq = 0.0;
    for (int i = 0; i < N; ++i) {
        const double d = noisy[i] - clean[i];
        noiseMax = std::max(noiseMax, std::fabs(d));
        noiseSq += d * d;
    }
    const double noiseRms = std::sqrt(noiseSq / N);
    std::cout << "MADDE 2 - GURULTU\n";
    std::cout << "  Gurultu orani = " << ratio << ", genlik = " << amp << " -> en buyuk gurultu <= " << ratio * amp << " (olculen en buyuk = " << noiseMax << ")\n";
    std::cout << "  Gurultu RMS: olculen = " << noiseRms << ", kuram (duzgun dagilim) = oran*genlik/sqrt(3) = " << ratio * amp / std::sqrt(3.0) << "\n";
    std::cout << "  Gurultulu sinyalin SNR'i (temiz sinyale gore) = " << std::setprecision(2) << snrDb(clean, noisy) << " dB\n\n" << std::setprecision(4);
    {
        const Series c{"temiz", {}, clean, green, true, false, false, 2};
        const Series n{"gurultulu", {}, noisy, red, true, true, false, 1};
        const double yN = amp * (1.0 + ratio) + 0.2;
        cv::Mat p = sinyal::plotSeries({c, n}, "Madde 2: gurultulu sinus (gurultu orani = " + sinyal::fmt(ratio, 2) + ")", -yN, yN, -1, N, "eleman no (n)");
        cv::imwrite((outDir / "2_gurultulu_sinyal.png").string(), p);
    }

    // =========================================================================================
    // MADDE 3: Gauss ve medyan filtre
    // =========================================================================================
    const int fSize = 9;
    const double fSigma = 2.0;
    const std::vector<double> gauss = sinyal::gaussianFilter(noisy, fSize, fSigma);
    const std::vector<double> med = sinyal::medianFilter(noisy, fSize);

    std::cout << "MADDE 3 - FILTRELEME (pencere = " << fSize << ", Gauss sigma = " << fSigma << ")\n";
    std::cout << "  " << std::left << std::setw(20) << "Sinyal" << std::right << std::setw(12) << "MSE" << std::setw(12) << "SNR (dB)" << std::setw(26) << "genlik (" + sinyal::fmt(freq, 0) + " Hz bileseni)" << "\n";
    auto row = [&](const std::string& name, const std::vector<double>& v) {
        std::cout << "  " << std::left << std::setw(20) << name << std::right << std::setw(12) << std::setprecision(4) << mse(v, clean) << std::setw(12)
                  << std::setprecision(2) << snrDb(clean, v) << std::setw(26) << std::setprecision(4) << ampAtBin(v, bin) << "\n";
    };
    row("gurultulu", noisy);
    row("Gauss sonrasi", gauss);
    row("Medyan sonrasi", med);
    std::cout << "  (temiz sinyalin " << freq << " Hz bilesen genligi = " << ampAtBin(clean, bin) << ")\n\n";

    {
        const Series c{"temiz", {}, clean, green, true, false, false, 2};
        const Series n{"gurultulu", {}, noisy, gray, true, false, false, 1};
        const Series g{"Gauss sonrasi", {}, gauss, blue, true, false, false, 2};
        const Series m{"Medyan sonrasi", {}, med, red, true, false, false, 2};
        const double yN = amp * (1.0 + ratio) + 0.2;
        cv::Mat p1 = sinyal::plotSeries({n, c, g}, "Madde 3: Gauss filtresi (pencere " + std::to_string(fSize) + ", sigma " + sinyal::fmt(fSigma, 1) + ")", -yN, yN, -1, N, "eleman no (n)");
        cv::Mat p2 = sinyal::plotSeries({n, c, m}, "Madde 3: Medyan filtresi (pencere " + std::to_string(fSize) + ")", -yN, yN, -1, N, "eleman no (n)");
        cv::Mat both;
        cv::vconcat(std::vector<cv::Mat>{p1, p2}, both);
        cv::imwrite((outDir / "3_gauss_ve_medyan.png").string(), both);
    }

    // Gurultu orani degisince filtreler
    std::cout << "  Gurultu oranina gore MSE (pencere " << fSize << "):\n";
    std::cout << "    " << std::setw(8) << "oran" << std::setw(14) << "gurultulu" << std::setw(14) << "Gauss" << std::setw(14) << "Medyan" << "\n";
    for (double r : {0.1, 0.3, 0.6, 1.0}) {
        const std::vector<double> nz = sinyal::addNoise(clean, r, 42);
        std::cout << "    " << std::setw(8) << std::setprecision(1) << r << std::setprecision(4) << std::setw(14) << mse(nz, clean) << std::setw(14)
                  << mse(sinyal::gaussianFilter(nz, fSize, fSigma), clean) << std::setw(14) << mse(sinyal::medianFilter(nz, fSize), clean) << "\n";
    }

    // Pencere boyutuna gore (gurultu orani = ratio)
    std::cout << "\n  Pencere boyutuna gore (gurultu orani " << std::setprecision(2) << ratio << ", Gauss sigma = pencere/4): MSE ve " << freq << " Hz bilesen genligi\n";
    std::cout << "    " << std::setw(8) << "pencere" << std::setw(14) << "Gauss MSE" << std::setw(14) << "Gauss gen." << std::setw(14) << "Medyan MSE" << std::setw(14) << "Medyan gen." << "\n";
    for (int w : {3, 5, 9, 15, 25, 41}) {
        const std::vector<double> g = sinyal::gaussianFilter(noisy, w, w / 4.0);
        const std::vector<double> m = sinyal::medianFilter(noisy, w);
        std::cout << "    " << std::setw(8) << w << std::setprecision(4) << std::setw(14) << mse(g, clean) << std::setw(14) << ampAtBin(g, bin) << std::setw(14) << mse(m, clean)
                  << std::setw(14) << ampAtBin(m, bin) << "\n";
    }

    // Ek deney: durtu (impuls) gurultusu
    std::cout << "\n  EK DENEY - durtu gurultusu (elemanlarin bir kismi +-1.5*genlik'e ayarlanir), pencere " << fSize << ":\n";
    std::cout << "    " << std::setw(10) << "oran" << std::setw(14) << "gurultulu" << std::setw(14) << "Gauss" << std::setw(14) << "Medyan" << "\n";
    for (double f : {0.02, 0.05, 0.10}) {
        const std::vector<double> imp = sinyal::addImpulseNoise(clean, f, 7);
        std::cout << "    " << std::setw(10) << std::setprecision(2) << f << std::setprecision(4) << std::setw(14) << mse(imp, clean) << std::setw(14)
                  << mse(sinyal::gaussianFilter(imp, fSize, fSigma), clean) << std::setw(14) << mse(sinyal::medianFilter(imp, fSize), clean) << "\n";
    }
    std::cout << "\n";

    // =========================================================================================
    // MADDE 4: DFT
    // =========================================================================================
    std::cout << "MADDE 4 - DFT (harmonik vektoru: k = 0.." << N / 2 << ", genlik cinsinden; k. harmonik = " << fs / N << " Hz x k)\n";
    // dogrulama: ters DFT temiz vektoru geri vermeli
    {
        const std::vector<double> back = sinyal::idftReal(sinyal::dftComplex(noisy));
        double maxErr = 0.0;
        for (int i = 0; i < N; ++i) maxErr = std::max(maxErr, std::fabs(back[i] - noisy[i]));
        std::cout << "  Dogrulama: ters DFT ile geri elde edilen vektor ile orijinal arasindaki en buyuk fark = " << std::scientific << maxErr << std::fixed << "\n";
    }

    struct Item {
        std::string name;
        const std::vector<double>* v;
    };
    const std::vector<Item> items = {{"temiz", &clean}, {"gurultulu", &noisy}, {"Gauss sonrasi", &gauss}, {"Medyan sonrasi", &med}};
    std::vector<std::vector<double>> amps(items.size());
    std::cout << "  " << std::left << std::setw(16) << "Sinyal" << std::right << std::setw(10) << "tepe k" << std::setw(14) << "tepe genlik" << std::setw(14) << "faz (derece)" << std::setw(16)
              << "taban ortalama" << std::setw(14) << "taban en buyuk" << std::setw(14) << "tepe/taban" << "\n";
    for (size_t i = 0; i < items.size(); ++i) {
        const SpecInfo s = analyse(*items[i].v, &amps[i]);
        std::cout << "  " << std::left << std::setw(16) << items[i].name << std::right << std::setw(10) << s.peakBin << std::setprecision(4) << std::setw(14) << s.peakAmp
                  << std::setprecision(2) << std::setw(14) << s.phaseDeg << std::setprecision(5) << std::setw(16) << s.floorMean << std::setw(14) << s.floorMax
                  << std::setprecision(1) << std::setw(14) << (s.floorMean > 0 ? s.peakAmp / s.floorMean : 0.0) << "\n";
    }
    std::cout << "  (gercek: tepe k = " << bin << ", genlik = " << std::setprecision(4) << amp << ", faz = " << std::setprecision(2) << phaseDeg << " derece)\n";
    {
        // gurultulu vektorun en guclu 5 harmonigi
        std::vector<std::pair<double, int>> top;
        for (size_t k = 0; k < amps[1].size(); ++k) top.push_back({amps[1][k], static_cast<int>(k)});
        std::sort(top.begin(), top.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
        std::cout << "  Gurultulu vektorun en guclu 5 harmonigi:";
        for (int i = 0; i < 5 && i < static_cast<int>(top.size()); ++i) std::cout << " k=" << top[i].second << " (" << std::setprecision(4) << top[i].first << ")";
        std::cout << "\n\n";
    }

    {
        const double yS = amp * 1.2;
        std::vector<double> ks(amps[0].size());
        for (size_t k = 0; k < ks.size(); ++k) ks[k] = static_cast<double>(k) * fs / N;
        const cv::Scalar cols[4] = {green, red, blue, orange};
        std::vector<cv::Mat> panels;
        for (size_t i = 0; i < items.size(); ++i) {
            const Series s{items[i].name + " (harmonik genlikleri)", ks, amps[i], cols[i], false, true, true, 1};
            panels.push_back(sinyal::plotSeries({s}, "Madde 4: DFT harmonik vektoru - " + items[i].name, 0.0, yS, 0, fs / 2.0, "frekans (Hz)", 1100, 300));
        }
        cv::Mat all;
        cv::vconcat(panels, all);
        cv::imwrite((outDir / "4_dft_harmonikler.png").string(), all);
    }

    std::cout << "Kaydedilenler (" << outDir.string() << "):\n"
              << "  1_sinus_vektoru.png     : Madde 1, vektor elemanlari (noktalar) ve surekli sinus dalgasi\n"
              << "  2_gurultulu_sinyal.png  : Madde 2, temiz ve gurultulu sinyal\n"
              << "  3_gauss_ve_medyan.png   : Madde 3, ust: Gauss, alt: Medyan sonucu\n"
              << "  4_dft_harmonikler.png   : Madde 4, DFT harmonik genlikleri (yukaridan asagi: temiz, gurultulu, Gauss, Medyan)\n";
    return 0;
}
