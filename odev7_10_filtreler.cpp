// Ödev 7 - Görev 10: Frekans uzayında low-pass (alçak geçiren), high-pass (yüksek geçiren), band-pass (bant geçiren)
// ve band-reject (bant durduran) filtreleri; bir resim ve gürültülü halleri üzerinde kullanılır; gürültüye, kenar
// geçişlerine ve genel parlaklığa etkisi ölçülür.
//
// Filtreler, resmin 2 boyutlu DFT'sindeki her harmonik r = sqrt(fu^2 + fv^2) (sıfır frekanstan uzaklık) yarıçapına göre
// bir kazançla (maske) çarpılarak uygulanır, sonra ters DFT yapılır.
//   Ideal:
//     LP  : r <= R ise 1, değilse 0                  HP : r >= R ise 1, değilse 0
//     BP  : R1 <= r <= R2 ise 1, değilse 0           BR : 1 - BP
//   Gauss (yumuşak geçişli):
//     LP  : exp(-r^2 / (2 s^2)), s = R / 1.1774 (r = R'de kazanç 0.5)      HP : 1 - LP
//     BP  : exp(-((r^2 - r0^2) / (r * Wb))^2), r0 = (R1+R2)/2, Wb = R2 - R1 (r = 0'da 0)      BR : 1 - BP
// Sıfır frekans (DC) resmin ortalama parlaklığıdır: LP ve BR bunu korur, HP ve BP sıfırlar.
//
// Resim: verilen yolda varsa o kullanılır; yoksa halka desenli sentetik bir "zone plate" üretilir. 256x256'ya getirilir.
//
// Kullanım: ./filtreler [cikti_dizini] [resim_yolu=images/zoneplate.png] [R=24] [R1=20] [R2=50] [sigma=20] [tuz_biber_orani=0.05]
// Derleme:  g++ -O2 odev7_10_filtreler.cpp -o filtreler.exe -std=c++17 $(pkg-config --cflags --libs opencv5)
// (sinyal.hpp, dft2d.hpp, dft2d_filtre.hpp, gurultu_resim.hpp ve metrikler.hpp bu dosyayla aynı klasörde olmalıdır)

#include "dft2d_filtre.hpp"
#include "gurultu_resim.hpp"
#include "metrikler.hpp"

#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fsys = std::filesystem;

// ---------------------------------------------------------------------------------------------
// Maskeler ve filtreleme
// ---------------------------------------------------------------------------------------------
static cv::Mat makeMask(int H, int W, const std::function<double(double)>& gain) {
    cv::Mat m(H, W, CV_64F);
    for (int v = 0; v < H; ++v)
        for (int u = 0; u < W; ++u) m.at<double>(v, u) = gain(sinyal::freqRadius(u, v, W, H));
    return m;
}

static cv::Mat applyMask(const sinyal::Spectrum2D& S, const cv::Mat& mask) {
    sinyal::Spectrum2D Y = S;
    for (int v = 0; v < S.H; ++v)
        for (int u = 0; u < S.W; ++u) Y.X[static_cast<size_t>(v) * S.W + u] *= mask.at<double>(v, u);
    return sinyal::idft2dReal(Y);  // kirpilmamis gercek degerler (CV_64F)
}

// ---------------------------------------------------------------------------------------------
// Olcumler (CV_64F resimler)
// ---------------------------------------------------------------------------------------------
static double meanOf(const cv::Mat& f) { return cv::mean(f)[0]; }

static double stdOf(const cv::Mat& f) {
    cv::Scalar m, s;
    cv::meanStdDev(f, m, s);
    return s[0];
}

static double edgeStrength(const cv::Mat& f) {  // ortalama |yatay fark| + |dikey fark|
    double s = 0.0;
    for (int y = 0; y + 1 < f.rows; ++y)
        for (int x = 0; x + 1 < f.cols; ++x) s += std::fabs(f.at<double>(y, x + 1) - f.at<double>(y, x)) + std::fabs(f.at<double>(y + 1, x) - f.at<double>(y, x));
    return s / (static_cast<double>(f.rows - 1) * (f.cols - 1));
}

static double rmsOf(const cv::Mat& f) {
    double s = 0.0;
    for (int y = 0; y < f.rows; ++y)
        for (int x = 0; x < f.cols; ++x) s += f.at<double>(y, x) * f.at<double>(y, x);
    return std::sqrt(s / (static_cast<double>(f.rows) * f.cols));
}

static double meanSquareMask(const cv::Mat& mask) {
    double s = 0.0;
    for (int v = 0; v < mask.rows; ++v)
        for (int u = 0; u < mask.cols; ++u) s += mask.at<double>(v, u) * mask.at<double>(v, u);
    return s / (static_cast<double>(mask.rows) * mask.cols);
}

static std::string f(double v, int p) { return sinyal::fmt(v, p); }

// Gosterim: HP ve BP ciktilari sifir ortalamali oldugu icin +128 ofsetle gosterilir
static cv::Mat showImage(const cv::Mat& f64, bool offset128) {
    cv::Mat out(f64.size(), CV_8UC1);
    for (int y = 0; y < f64.rows; ++y)
        for (int x = 0; x < f64.cols; ++x) out.at<uchar>(y, x) = cv::saturate_cast<uchar>(f64.at<double>(y, x) + (offset128 ? 128.0 : 0.0));
    return out;
}

static cv::Mat labelled(const cv::Mat& g, const std::string& l1, const std::string& l2 = "") {
    cv::Mat bgr;
    cv::cvtColor(g, bgr, cv::COLOR_GRAY2BGR);
    cv::putText(bgr, l1, cv::Point(5, 16), cv::FONT_HERSHEY_SIMPLEX, 0.42, cv::Scalar(0, 255, 255), 1, cv::LINE_AA);
    if (!l2.empty()) cv::putText(bgr, l2, cv::Point(5, 33), cv::FONT_HERSHEY_SIMPLEX, 0.42, cv::Scalar(0, 255, 255), 1, cv::LINE_AA);
    return bgr;
}

static cv::Mat makeZonePlate(int N) {
    cv::Mat z = cv::Mat::zeros(N, N, CV_8UC1);
    const double c = (N - 1) / 2.0, rMax = 0.97 * N / 2.0, s = 248.0 * N / 256.0;  // s: yerel frekans r/s dongu/piksel
    for (int y = 0; y < N; ++y)
        for (int x = 0; x < N; ++x) {
            const double r2 = (x - c) * (x - c) + (y - c) * (y - c);
            if (r2 <= rMax * rMax) z.at<uchar>(y, x) = cv::saturate_cast<uchar>(255.0 * 0.5 * (1.0 + std::cos(sinyal::PI * r2 / s)));
        }
    return z;
}

struct NamedFilter {
    std::string name;
    std::string kind;  // LP, HP, BP, BR
    bool ideal;
    cv::Mat mask;
};

int main(int argc, char** argv) {
    const fsys::path outDir = (argc >= 2) ? fsys::path(argv[1]) : fsys::path("cikti_filtreler");
    const std::string path = (argc >= 3) ? argv[2] : "images/zoneplate.png";
    const double R = (argc >= 4) ? std::stod(argv[3]) : 24.0;
    const double R1 = (argc >= 5) ? std::stod(argv[4]) : 20.0;
    const double R2 = (argc >= 6) ? std::stod(argv[5]) : 50.0;
    const double sigma = (argc >= 7) ? std::stod(argv[6]) : 20.0;
    const double ratio = (argc >= 8) ? std::stod(argv[7]) : 0.05;
    fsys::create_directories(outDir);
    const int N = 256;

    cv::Mat orig;
    {
        cv::Mat src = cv::imread(path, cv::IMREAD_GRAYSCALE);
        if (src.empty()) {
            std::cout << "(" << path << " bulunamadi: sentetik zone plate uretildi)\n";
            orig = makeZonePlate(N);
            cv::imwrite((outDir / "10_sentetik_zoneplate.png").string(), orig);
        } else {
            cv::resize(src, orig, cv::Size(N, N), 0, 0, cv::INTER_AREA);
            std::cout << "Resim: " << path << " (" << src.cols << "x" << src.rows << " -> " << N << "x" << N << ")\n";
        }
    }
    const cv::Mat gauss = sinyal::addGaussianNoise(orig, sigma, 1);
    const cv::Mat sp = sinyal::addSaltAndPepperNoise(orig, ratio, 2);

    // ---- filtreler ----
    const double sG = R / 1.1774;  // Gauss LP: r = R'de kazanc 0.5
    const double r0 = 0.5 * (R1 + R2), Wb = R2 - R1;
    auto bpGauss = [&](double r) { return r < 1e-9 ? 0.0 : std::exp(-std::pow((r * r - r0 * r0) / (r * Wb), 2.0)); };
    auto lpGauss = [&](double r) { return std::exp(-r * r / (2.0 * sG * sG)); };
    std::vector<NamedFilter> filters;
    filters.push_back({"Ideal LP  R=" + f(R, 0), "LP", true, makeMask(N, N, [&](double r) { return r <= R ? 1.0 : 0.0; })});
    filters.push_back({"Ideal HP  R=" + f(R, 0), "HP", true, makeMask(N, N, [&](double r) { return r >= R ? 1.0 : 0.0; })});
    filters.push_back({"Ideal BP  " + f(R1, 0) + "-" + f(R2, 0), "BP", true, makeMask(N, N, [&](double r) { return (r >= R1 && r <= R2) ? 1.0 : 0.0; })});
    filters.push_back({"Ideal BR  " + f(R1, 0) + "-" + f(R2, 0), "BR", true, makeMask(N, N, [&](double r) { return (r >= R1 && r <= R2) ? 0.0 : 1.0; })});
    filters.push_back({"Gauss LP  R=" + f(R, 0), "LP", false, makeMask(N, N, lpGauss)});
    filters.push_back({"Gauss HP  R=" + f(R, 0), "HP", false, makeMask(N, N, [&](double r) { return 1.0 - lpGauss(r); })});
    filters.push_back({"Gauss BP  " + f(R1, 0) + "-" + f(R2, 0), "BP", false, makeMask(N, N, bpGauss)});
    filters.push_back({"Gauss BR  " + f(R1, 0) + "-" + f(R2, 0), "BR", false, makeMask(N, N, [&](double r) { return 1.0 - bpGauss(r); })});

    std::cout << "FILTRELER (r = sifir frekanstan uzaklik, 'dongu/resim'; Nyquist siniri = 128): LP/HP kesim R = " << R << ", BP/BR bandi " << R1 << " - " << R2 << "\n";
    std::cout << "  Gurultu: Gauss sigma = " << sigma << ", tuz-biber orani = " << ratio << "\n\n";

    // ---- spektrumlar ve filtrelenmis resimler ----
    const cv::Mat inputs[3] = {orig, gauss, sp};
    const std::string inNames[3] = {"TEMIZ RESIM", "GAUSS GURULTULU RESIM", "TUZ-BIBER GURULTULU RESIM"};
    sinyal::Spectrum2D spec[3];
    for (int i = 0; i < 3; ++i) spec[i] = sinyal::dft2dComplex(inputs[i]);
    std::vector<std::vector<cv::Mat>> out(3);  // out[giris][filtre], kirpilmamis
    for (int i = 0; i < 3; ++i)
        for (const NamedFilter& nf : filters) out[i].push_back(applyMask(spec[i], nf.mask));

    cv::Mat o64;
    orig.convertTo(o64, CV_64F);
    const double origStd = stdOf(o64), origEdge = edgeStrength(o64), origMean = meanOf(o64);

    // ===== Tablo 1: temiz resim =====
    std::cout << "TABLO 1 - TEMIZ RESIM UZERINDE FILTRELERIN ETKISI (orijinal: ort. parlaklik " << f(origMean, 2) << ", std " << f(origStd, 2) << ", kenar gucu " << f(origEdge, 3) << ")\n";
    std::cout << "  " << std::left << std::setw(20) << "Filtre" << std::right << std::setw(12) << "DC kazanci" << std::setw(13) << "ort. parlak." << std::setw(10) << "min" << std::setw(10) << "max" << std::setw(12) << "kenar gucu"
              << std::setw(10) << "oran" << std::setw(10) << "ZNCC" << std::setw(14) << "enerji orani" << "\n";
    for (size_t k = 0; k < filters.size(); ++k) {
        const cv::Mat& y = out[0][k];
        double mn, mx;
        cv::minMaxLoc(y, &mn, &mx);
        const double eStr = edgeStrength(y);
        // ZNCC: kirpilmis 8-bit orijinalle (HP/BP icin anlamli degil ama yine de verilir)
        const double zn = sinyal::znccImg(showImage(y, false), orig);
        std::cout << "  " << std::left << std::setw(20) << filters[k].name << std::right << std::setw(12) << f(filters[k].mask.at<double>(0, 0), 2) << std::setw(13) << f(meanOf(y), 2) << std::setw(10) << f(mn, 1)
                  << std::setw(10) << f(mx, 1) << std::setw(12) << f(eStr, 3) << std::setw(10) << f(eStr / origEdge, 2) << std::setw(10) << (std::isnan(zn) ? "-" : f(zn, 3)) << std::setw(14)
                  << f(stdOf(y) * stdOf(y) / (origStd * origStd), 3) << "\n";
    }
    std::cout << "  (ort. parlak. = filtre ciktisinin kirpilmamis ortalamasi = DC kazanci x orijinal ortalama; min/max kirpilmamis degerler: 0-255 disina tasma = halka (ringing) / Gibbs etkisi)\n"
              << "  (oran = filtre cikti kenar gucu / orijinal kenar gucu; enerji orani = cikti varyansi / orijinal varyans)\n\n";

    // ===== Tablo 2: gurultu =====
    for (int i = 1; i < 3; ++i) {
        cv::Mat noise = [&] {
            cv::Mat a, b;
            inputs[i].convertTo(a, CV_64F);
            return a - o64;
        }();
        const double nRms = rmsOf(noise);
        std::cout << "TABLO 2" << (i == 1 ? "a" : "b") << " - " << inNames[i] << ": filtrelerin gurultuye etkisi (giris gurultu RMS = " << f(nRms, 2) << ", giris SNR = " << f(20.0 * std::log10(origStd / nRms), 2) << " dB)\n";
        std::cout << "  " << std::left << std::setw(20) << "Filtre" << std::right << std::setw(18) << "gecen gurultu RMS" << std::setw(16) << "kuram (beyaz)" << std::setw(14) << "gecen oran" << std::setw(14) << "cikis SNR dB" << std::setw(14) << "SNR degisimi" << "\n";
        for (size_t k = 0; k < filters.size(); ++k) {
            cv::Mat passed = out[i][k] - out[0][k];  // filtre(gurultulu) - filtre(temiz) = filtre(gurultu) (dogrusal)
            const double pRms = rmsOf(passed);
            const double theory = nRms * std::sqrt(meanSquareMask(filters[k].mask));
            const double sigOut = stdOf(out[0][k]);
            const double snrOut = sigOut > 1e-9 && pRms > 1e-9 ? 20.0 * std::log10(sigOut / pRms) : 999.0;
            const double snrIn = 20.0 * std::log10(origStd / nRms);
            std::cout << "  " << std::left << std::setw(20) << filters[k].name << std::right << std::setw(18) << f(pRms, 2) << std::setw(16) << f(theory, 2) << std::setw(14) << f(pRms / nRms, 3) << std::setw(14)
                      << f(snrOut, 2) << std::setw(14) << f(snrOut - snrIn, 2) << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "  (gecen gurultu = filtre(gurultulu resim) - filtre(temiz resim); kuram = giris gurultu RMS x sqrt(maskenin karelerinin ortalamasi), beyaz gurultu varsayimi;\n"
              << "   cikis SNR = filtrelenmis temiz resmin std'si / gecen gurultu RMS, dB)\n\n";

    // ---- gorseller ----
    auto montageFor = [&](int i) {
        cv::Mat rowIdeal, rowGauss;
        std::vector<cv::Mat> pi = {labelled(inputs[i], i == 0 ? "temiz resim" : (i == 1 ? "Gauss gurultulu" : "tuz-biber gurultulu"))};
        std::vector<cv::Mat> pg = {labelled(inputs[i], i == 0 ? "temiz resim" : (i == 1 ? "Gauss gurultulu" : "tuz-biber gurultulu"))};
        for (size_t k = 0; k < filters.size(); ++k) {
            const bool off = (filters[k].kind == "HP" || filters[k].kind == "BP");
            cv::Mat p = labelled(showImage(out[i][k], off), filters[k].name, off ? "(+128 ofsetli gosterim)" : "");
            (filters[k].ideal ? pi : pg).push_back(p);
        }
        cv::hconcat(pi, rowIdeal);
        cv::hconcat(pg, rowGauss);
        cv::Mat m;
        cv::vconcat(std::vector<cv::Mat>{rowIdeal, rowGauss}, m);
        return m;
    };
    cv::imwrite((outDir / "10_temiz_filtreler.png").string(), montageFor(0));
    cv::imwrite((outDir / "10_gauss_gurultulu_filtreler.png").string(), montageFor(1));
    cv::imwrite((outDir / "10_tuz_biber_filtreler.png").string(), montageFor(2));

    // maskeler (sifir frekans ortada) ve orijinalin spektrumu
    {
        std::vector<cv::Mat> top, bottom;
        for (size_t k = 0; k < filters.size(); ++k) {
            cv::Mat m8(N, N, CV_8UC1);
            for (int v = 0; v < N; ++v)
                for (int u = 0; u < N; ++u) m8.at<uchar>((v + N / 2) % N, (u + N / 2) % N) = cv::saturate_cast<uchar>(255.0 * filters[k].mask.at<double>(v, u));
            (filters[k].ideal ? top : bottom).push_back(labelled(m8, filters[k].name, "filtre maskesi (DC ortada)"));
        }
        cv::Mat r1, r2, all;
        cv::hconcat(top, r1);
        cv::hconcat(bottom, r2);
        cv::vconcat(std::vector<cv::Mat>{r1, r2}, all);
        cv::imwrite((outDir / "10_filtre_maskeleri.png").string(), all);

        const cv::Mat amp = sinyal::harmonicsFromSpectrum(spec[0]);
        cv::imwrite((outDir / "10_orijinal_spektrum.png").string(), labelled(sinyal::spectrumToImage(amp), "orijinal resmin log spektrumu", "(DC ortada)"));
    }

    std::cout << "Kaydedilenler (" << outDir.string() << "):\n"
              << "  10_temiz_filtreler.png, 10_gauss_gurultulu_filtreler.png, 10_tuz_biber_filtreler.png : ust satir ideal, alt satir Gauss filtreler; soldan saga girdi, LP, HP, BP, BR\n"
              << "  10_filtre_maskeleri.png                                                      : filtrelerin frekans maskeleri\n"
              << "  10_orijinal_spektrum.png                                                     : orijinal resmin log spektrumu\n";
    return 0;
}
