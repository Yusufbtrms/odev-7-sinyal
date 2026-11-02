// Ödev 7 - Görev 5 ve 6: Birden çok sinüsün toplamı, DFT analizi ve frekans filtresi.
//
// Görev 5: Frekansı 50 Hz genliği 12 olan ve frekansı 30 Hz genliği 5 olan, gürültü eklenmiş iki sinyalin önce
//          ayrı ayrı DFT analizi yapılır. Sonra iki sinyal toplanır; toplam sinyal ve harmonik vektörü çizilir.
// Görev 6: Harmonik vektöründe eşiğin altında kalan harmonikler atılarak ve ters DFT ile orijinal sinyal yeniden
//          üretilir (frekans filtresi). Görev 5'teki iki sinyal ve toplamları için test edilir.
//
// Sinyal: N = 256 örnek, örnekleme frekansı 256 Hz (vektör 1 saniye, her harmonik 1 Hz; 30 ve 50 Hz tam harmoniklerdir).
//
// Kullanım: ./frekans [cikti_dizini] [gurultu_orani=0.3] [esik=1.5] [N=256]
// Derleme:  g++ -O2 odev7_5_6_frekans.cpp -o frekans.exe -std=c++17 $(pkg-config --cflags --libs opencv5)
// (sinyal.hpp ve frekans.hpp bu dosyayla aynı klasörde olmalıdır)

#include "frekans.hpp"

#include <opencv2/imgcodecs.hpp>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fsys = std::filesystem;
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

static double rmsDiff(const std::vector<double>& a, const std::vector<double>& b) { return std::sqrt(mse(a, b)); }

// sin(w n + faz) icin DFT fazi = faz - 90 derece; bu yuzden tahmin = arg + 90 derece
static double phaseEstimateDeg(const std::vector<double>& phase, int bin) {
    double ph = phase[bin] + sinyal::PI / 2.0;
    while (ph > sinyal::PI) ph -= 2.0 * sinyal::PI;
    while (ph < -sinyal::PI) ph += 2.0 * sinyal::PI;
    return ph * 180.0 / sinyal::PI;
}

static double floorMean(const std::vector<double>& amp, const std::vector<int>& excluded, double* floorMaxOut = nullptr) {
    double sum = 0.0, mx = 0.0;
    int cnt = 0;
    for (size_t k = 1; k < amp.size(); ++k) {
        if (std::find(excluded.begin(), excluded.end(), static_cast<int>(k)) != excluded.end()) continue;
        sum += amp[k];
        mx = std::max(mx, amp[k]);
        ++cnt;
    }
    if (floorMaxOut) *floorMaxOut = mx;
    return cnt ? sum / cnt : 0.0;
}

int main(int argc, char** argv) {
    const fsys::path outDir = (argc >= 2) ? fsys::path(argv[1]) : fsys::path("cikti_frekans");
    const double ratio = (argc >= 3) ? std::stod(argv[2]) : 0.3;
    const double threshold = (argc >= 4) ? std::stod(argv[3]) : 1.5;
    const int N = (argc >= 5) ? std::stoi(argv[4]) : 256;
    fsys::create_directories(outDir);
    if (!sinyal::isPowerOfTwo(static_cast<size_t>(N))) {
        std::cerr << "Hata: N 2'nin kuvveti olmali\n";
        return 1;
    }
    const double fsHz = static_cast<double>(N);
    const double f1 = 50.0, a1 = 12.0, f2 = 30.0, a2 = 5.0;
    const int bin1 = static_cast<int>(std::lround(f1 * N / fsHz)), bin2 = static_cast<int>(std::lround(f2 * N / fsHz));

    // ---- sinyaller ----
    std::vector<double> c1(N), c2(N);
    sinyal::fillSine(c1, a1, 0.0, f1, fsHz);
    sinyal::fillSine(c2, a2, 0.0, f2, fsHz);
    const std::vector<double> n1 = sinyal::addNoise(c1, ratio, 42);  // sinyal 1 + gurultu
    const std::vector<double> n2 = sinyal::addNoise(c2, ratio, 43);  // sinyal 2 + gurultu (bagimsiz gurultu)
    const std::vector<double> cSum = sinyal::addSignals(c1, c2);
    const std::vector<double> nSum = sinyal::addSignals(n1, n2);

    const cv::Scalar blue(200, 120, 0), red(0, 0, 220), green(0, 150, 0), gray(170, 170, 170), orange(0, 140, 255);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "SINYALLER: N = " << N << ", ornekleme = " << fsHz << " Hz, gurultu orani = " << std::setprecision(2) << ratio << std::setprecision(4) << "\n";
    std::cout << "  Sinyal 1: " << f1 << " Hz, genlik " << a1 << " (harmonik k = " << bin1 << ")\n";
    std::cout << "  Sinyal 2: " << f2 << " Hz, genlik " << a2 << " (harmonik k = " << bin2 << ")\n\n";

    // =========================================================================================
    // GOREV 5: ayri ve toplam DFT analizi
    // =========================================================================================
    std::cout << "GOREV 5 - DFT ANALIZI (gurultulu sinyaller)\n";
    std::vector<double> ph1, ph2, phS;
    const std::vector<double> amp1 = sinyal::dftHarmonics(n1, &ph1);
    const std::vector<double> amp2 = sinyal::dftHarmonics(n2, &ph2);
    const std::vector<double> ampS = sinyal::dftHarmonics(nSum, &phS);
    const std::vector<double> ampSClean = sinyal::dftHarmonics(cSum);

    auto peakOf = [](const std::vector<double>& amp) {
        int pk = 1;
        for (size_t k = 1; k < amp.size(); ++k)
            if (amp[k] > amp[pk]) pk = static_cast<int>(k);
        return pk;
    };
    std::cout << "  " << std::left << std::setw(26) << "Sinyal" << std::right << std::setw(9) << "tepe k" << std::setw(14) << "tepe genlik" << std::setw(14) << "faz (derece)"
              << std::setw(16) << "taban ortalama" << std::setw(15) << "taban en buyuk" << "\n";
    {
        const int p = peakOf(amp1);
        double fm = 0;
        const double fl = floorMean(amp1, {p}, &fm);
        std::cout << "  " << std::left << std::setw(26) << "sinyal 1 (50 Hz) gurultulu" << std::right << std::setw(9) << p << std::setw(14) << amp1[p] << std::setprecision(2)
                  << std::setw(14) << phaseEstimateDeg(ph1, p) << std::setprecision(4) << std::setw(16) << fl << std::setw(15) << fm << "\n";
    }
    {
        const int p = peakOf(amp2);
        double fm = 0;
        const double fl = floorMean(amp2, {p}, &fm);
        std::cout << "  " << std::left << std::setw(26) << "sinyal 2 (30 Hz) gurultulu" << std::right << std::setw(9) << p << std::setw(14) << amp2[p] << std::setprecision(2)
                  << std::setw(14) << phaseEstimateDeg(ph2, p) << std::setprecision(4) << std::setw(16) << fl << std::setw(15) << fm << "\n";
    }
    std::cout << "\n  TOPLAM sinyal (gurultulu sinyal 1 + gurultulu sinyal 2):\n";
    {
        std::vector<std::pair<double, int>> top;
        for (size_t k = 1; k < ampS.size(); ++k) top.push_back({ampS[k], static_cast<int>(k)});
        std::sort(top.begin(), top.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
        std::cout << "    En guclu 5 harmonik:";
        for (int i = 0; i < 5; ++i) std::cout << " k=" << top[i].second << " (" << top[i].first << ")";
        std::cout << "\n";
        double fm = 0;
        const double fl = floorMean(ampS, {bin1, bin2}, &fm);
        std::cout << "    k=" << bin1 << ": genlik " << ampS[bin1] << ", faz " << std::setprecision(2) << phaseEstimateDeg(phS, bin1) << std::setprecision(4) << " derece; k=" << bin2 << ": genlik "
                  << ampS[bin2] << ", faz " << std::setprecision(2) << phaseEstimateDeg(phS, bin2) << std::setprecision(4) << " derece\n";
        const double sigmaN = rmsDiff(nSum, cSum);
        std::cout << "    Gurultu tabani (iki tepe disindaki harmoniklerin ortalamasi) = " << fl << ", en buyuk = " << fm << "; kuram: gurultu RMS * sqrt(pi/N) = " << sigmaN * std::sqrt(sinyal::PI / N)
                  << " (gurultu RMS = " << sigmaN << ")\n";
        std::cout << "    Toplam sinyalin zaman uzayinda SNR'i = " << std::setprecision(2) << snrDb(cSum, nSum) << " dB, frekans uzayinda zayif tepe (30 Hz) / taban ortalamasi = " << std::setprecision(1)
                  << ampS[bin2] / fl << " kat\n" << std::setprecision(4);
        // dogrusallik: DFT(toplam) = DFT(sinyal 1) + DFT(sinyal 2)
        const auto X1 = sinyal::dftComplex(n1), X2 = sinyal::dftComplex(n2), XS = sinyal::dftComplex(nSum);
        double maxDiff = 0.0;
        for (int k = 0; k < N; ++k) maxDiff = std::max(maxDiff, std::abs(XS[k] - X1[k] - X2[k]));
        std::cout << "    Dogrusallik kontrolu: |DFT(toplam) - DFT(sinyal 1) - DFT(sinyal 2)| en buyuk degeri = " << std::scientific << maxDiff << std::fixed << "\n\n";
    }

    // grafikler
    auto timePlot = [&](const std::string& title, const std::vector<Series>& s, double yR, double xMax) {
        return sinyal::plotSeries(s, title, -yR, yR, -1, xMax, "eleman no (n)");
    };
    auto freqAxis = [&](const std::vector<double>& amp) {
        std::vector<double> ks(amp.size());
        for (size_t k = 0; k < ks.size(); ++k) ks[k] = static_cast<double>(k) * fsHz / N;
        return ks;
    };
    auto specPlot = [&](const std::string& title, const std::vector<double>& amp, const cv::Scalar& col, double yMax) {
        const Series s{"harmonik genlikleri", freqAxis(amp), amp, col, false, true, true, 1};
        return sinyal::plotSeries({s}, title, 0.0, yMax, 0, fsHz / 2.0, "frekans (Hz)", 1100, 320);
    };
    {
        const double y1 = a1 * (1.0 + ratio) + 1.0;
        cv::Mat t1 = timePlot("Gorev 5: sinyal 1 (50 Hz, genlik 12) gurultulu, ilk 64 eleman", {Series{"temiz", {}, c1, green, true, false, false, 2}, Series{"gurultulu", {}, n1, red, true, true, false, 1}}, y1, 64);
        cv::Mat s1 = specPlot("Gorev 5: sinyal 1'in DFT harmonik vektoru", amp1, red, a1 * 1.2);
        cv::Mat o1;
        cv::vconcat(std::vector<cv::Mat>{t1, s1}, o1);
        cv::imwrite((outDir / "5_sinyal1_dft.png").string(), o1);

        const double y2 = a2 * (1.0 + ratio) + 1.0;
        cv::Mat t2 = timePlot("Gorev 5: sinyal 2 (30 Hz, genlik 5) gurultulu, ilk 64 eleman", {Series{"temiz", {}, c2, green, true, false, false, 2}, Series{"gurultulu", {}, n2, red, true, true, false, 1}}, y2, 64);
        cv::Mat s2 = specPlot("Gorev 5: sinyal 2'nin DFT harmonik vektoru", amp2, red, a1 * 1.2);
        cv::Mat o2;
        cv::vconcat(std::vector<cv::Mat>{t2, s2}, o2);
        cv::imwrite((outDir / "5_sinyal2_dft.png").string(), o2);

        const double yS = (a1 + a2) * (1.0 + ratio) + 1.0;
        const Series cS{"temiz toplam", {}, cSum, green, true, false, false, 2};
        const Series nS{"gurultulu toplam", {}, nSum, red, true, false, false, 1};
        cv::Mat tFull = timePlot("Gorev 5: toplam sinyal (tum vektor)", {nS, cS}, yS, N);
        cv::Mat tZoom = timePlot("Gorev 5: toplam sinyal (ilk 64 eleman)", {Series{"temiz toplam", {}, cSum, green, true, false, false, 2}, Series{"gurultulu toplam", {}, nSum, red, true, true, false, 1}}, yS, 64);
        cv::Mat sS = specPlot("Gorev 5: toplam sinyalin DFT harmonik vektoru (30 Hz ve 50 Hz tepeleri)", ampS, blue, a1 * 1.2);
        cv::Mat oS;
        cv::vconcat(std::vector<cv::Mat>{tFull, tZoom, sS}, oS);
        cv::imwrite((outDir / "5_toplam_sinyal_ve_dft.png").string(), oS);
    }

    // =========================================================================================
    // GOREV 6: frekans filtresi
    // =========================================================================================
    std::cout << "GOREV 6 - FREKANS FILTRESI (esik = " << std::setprecision(2) << threshold << " genlik birimi)\n" << std::setprecision(4);
    {
        // dogrulama: esik 0 ise sinyal aynen geri gelmeli
        const std::vector<double> back = sinyal::frequencyFilter(nSum, 0.0);
        double maxErr = 0.0;
        for (int i = 0; i < N; ++i) maxErr = std::max(maxErr, std::fabs(back[i] - nSum[i]));
        std::cout << "  Dogrulama: esik = 0 iken filtre cikisi ile girdi arasindaki en buyuk fark = " << std::scientific << maxErr << std::fixed << "\n";
    }

    struct Case {
        std::string name;
        const std::vector<double>* clean;
        const std::vector<double>* noisy;
        const std::vector<double>* amp;
        std::string file;
        double yAmp;
    };
    const std::vector<Case> cases = {{"sinyal 1 (50 Hz)", &c1, &n1, &amp1, "6_sinyal1_frekans_filtresi.png", a1},
                                     {"sinyal 2 (30 Hz)", &c2, &n2, &amp2, "6_sinyal2_frekans_filtresi.png", a2},
                                     {"toplam sinyal", &cSum, &nSum, &ampS, "6_toplam_frekans_filtresi.png", a1 + a2}};
    std::vector<std::vector<double>> filtered;
    std::cout << "\n  " << std::left << std::setw(20) << "Sinyal" << std::right << std::setw(14) << "gurultulu MSE" << std::setw(14) << "filtreli MSE" << std::setw(14) << "gurultulu SNR"
              << std::setw(14) << "filtreli SNR" << std::setw(12) << "korunan k" << "\n";
    for (const Case& c : cases) {
        std::vector<int> kept;
        std::vector<double> keptAmp;
        filtered.push_back(sinyal::frequencyFilter(*c.noisy, threshold, &kept, &keptAmp));
        const std::vector<double>& f = filtered.back();
        std::cout << "  " << std::left << std::setw(20) << c.name << std::right << std::setw(14) << mse(*c.noisy, *c.clean) << std::setw(14) << mse(f, *c.clean) << std::setprecision(2) << std::setw(11)
                  << snrDb(*c.clean, *c.noisy) << " dB" << std::setw(11) << snrDb(*c.clean, f) << " dB" << std::setprecision(4) << std::setw(8) << kept.size() << "\n";
        std::cout << "      korunan harmonikler:";
        for (size_t i = 0; i < kept.size(); ++i) std::cout << " k=" << kept[i] << " (" << std::setprecision(3) << keptAmp[i] << ")";
        std::cout << std::setprecision(4) << "\n";
    }

    // Esik taramasi (toplam sinyal)
    std::cout << "\n  ESIK TARAMASI (toplam sinyal, gurultu orani " << std::setprecision(2) << ratio << "): filtre ciktisinin temiz toplam sinyale gore hatasi\n" << std::setprecision(4);
    std::cout << "    " << std::setw(8) << "esik" << std::setw(14) << "korunan k" << std::setw(14) << "MSE" << std::setw(12) << "SNR (dB)" << "\n";
    for (double th : {0.0, 0.1, 0.3, 0.6, 1.0, 1.5, 3.0, 5.5, 13.0}) {
        std::vector<int> kept;
        const std::vector<double> f = sinyal::frequencyFilter(nSum, th, &kept);
        std::cout << "    " << std::setw(8) << std::setprecision(2) << th << std::setprecision(4) << std::setw(14) << kept.size() << std::setw(14) << mse(f, cSum) << std::setprecision(2) << std::setw(12)
                  << snrDb(cSum, f) << std::setprecision(4) << "\n";
    }

    // Gurultu orani degisince (esik sabit)
    std::cout << "\n  GURULTU ORANINA GORE (esik = " << std::setprecision(2) << threshold << ", toplam sinyal):\n" << std::setprecision(4);
    std::cout << "    " << std::setw(8) << "oran" << std::setw(14) << "gurultulu MSE" << std::setw(14) << "filtreli MSE" << std::setw(14) << "korunan k" << "\n";
    for (double r : {0.1, 0.3, 0.6, 1.0}) {
        const std::vector<double> ns = sinyal::addSignals(sinyal::addNoise(c1, r, 42), sinyal::addNoise(c2, r, 43));
        std::vector<int> kept;
        const std::vector<double> f = sinyal::frequencyFilter(ns, threshold, &kept);
        std::cout << "    " << std::setw(8) << std::setprecision(1) << r << std::setprecision(4) << std::setw(14) << mse(ns, cSum) << std::setw(14) << mse(f, cSum) << std::setw(14) << kept.size() << "\n";
    }

    // Zaman uzayi filtreleriyle karsilastirma
    {
        const std::vector<double> g = sinyal::gaussianFilter(nSum, 9, 2.0);
        const std::vector<double> m = sinyal::medianFilter(nSum, 9);
        std::cout << "\n  KARSILASTIRMA (toplam sinyal, MSE; zaman uzayi filtreleri Odev 7'nin onceki bolumundeki Gauss ve medyan, pencere 9):\n";
        std::cout << "    gurultulu = " << mse(nSum, cSum) << ", Gauss = " << mse(g, cSum) << ", Medyan = " << mse(m, cSum) << ", Frekans filtresi (esik " << std::setprecision(2) << threshold
                  << std::setprecision(4) << ") = " << mse(filtered[2], cSum) << "\n";
    }

    // grafikler: zaman uzayi (temiz, gurultulu, filtreli) ve esikli spektrum
    for (size_t i = 0; i < cases.size(); ++i) {
        const Case& c = cases[i];
        std::vector<int> kept;
        sinyal::frequencyFilter(*c.noisy, threshold, &kept);
        const double yR = c.yAmp * (1.0 + ratio) + 1.0;
        cv::Mat t = timePlot("Gorev 6: " + c.name + " - frekans filtresi sonucu (ilk 64 eleman)",
                             {Series{"gurultulu", {}, *c.noisy, gray, true, false, false, 1}, Series{"temiz", {}, *c.clean, green, true, false, false, 2}, Series{"filtreli", {}, filtered[i], red, true, false, false, 2}},
                             yR, 64);

        std::vector<double> kx, ky, rx, ry;
        for (size_t k = 0; k < c.amp->size(); ++k) {
            const bool isKept = std::find(kept.begin(), kept.end(), static_cast<int>(k)) != kept.end();
            (isKept ? kx : rx).push_back(static_cast<double>(k) * fsHz / N);
            (isKept ? ky : ry).push_back((*c.amp)[k]);
        }
        const Series removed{"atilan harmonikler (esigin altinda)", rx, ry, gray, false, true, true, 1};
        const Series keptS{"korunan harmonikler", kx, ky, green, false, true, true, 1};
        const Series thr{"esik = " + sinyal::fmt(threshold, 2), {0.0, fsHz / 2.0}, {threshold, threshold}, orange, true, false, false, 2};
        cv::Mat s = sinyal::plotSeries({removed, keptS, thr}, "Gorev 6: " + c.name + " - gurultulu sinyalin harmonik vektoru ve esik", 0.0, a1 * 1.2, 0, fsHz / 2.0, "frekans (Hz)", 1100, 340);
        cv::Mat o;
        cv::vconcat(std::vector<cv::Mat>{t, s}, o);
        cv::imwrite((outDir / c.file).string(), o);
    }

    std::cout << "\nKaydedilenler (" << outDir.string() << "):\n"
              << "  5_sinyal1_dft.png, 5_sinyal2_dft.png      : Gorev 5, ayri sinyaller ve harmonik vektorleri\n"
              << "  5_toplam_sinyal_ve_dft.png                : Gorev 5, toplam sinyal (tum vektor ve ilk 64 eleman) ve harmonik vektoru\n"
              << "  6_sinyal1_frekans_filtresi.png            : Gorev 6, sinyal 1 (gurultulu, temiz, filtreli) ve esikli spektrum\n"
              << "  6_sinyal2_frekans_filtresi.png            : Gorev 6, sinyal 2\n"
              << "  6_toplam_frekans_filtresi.png             : Gorev 6, toplam sinyal\n";
    return 0;
}
