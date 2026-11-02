# Ödev 7: Sinüs, Gürültü, Filtreleme ve DFT Raporu

Kodlar `sinyal.hpp` (fonksiyonlar) ve `odev7_sinyal.cpp` (ana program) dosyalarındadır. Çıktılar `cikti_sinyal/` klasöründe, ekran çıktısı `sinyal_sonuclari.txt` dosyasındadır.

## 1. Yöntem

- **Madde 1, `fillSine`:** Genlik, faz, frekans (Hz) ve örnekleme frekansı alır, vektörü `genlik * sin(2*pi*frekans*n/örnekleme + faz)` ile doldurur; vektör boyutu 2'nin kuvveti değilse hata verir (periyotla vermek için `fillSineByPeriod` da vardır). Kullanılan sinyal: N = 256, örnekleme 256 Hz (vektör 1 saniye), 8 Hz, genlik 2, faz 30°. Vektör, vektör elemanlarını nokta, sürekli sinüs dalgasını çizgi olarak gösteren grafikle çizildi.
- **Madde 2, `addNoise`:** Vektörün genliğini (en büyük mutlak değer) A bulur, her elemana `oran * A * u` ekler; u, [-1, 1] aralığında düzgün dağılımlı rastgele sayıdır. Böylece en büyük gürültü oran × A kadardır. Oran [0, 1] dışındaysa hata verir. Kullanılan oran: 0.3.
- **Madde 3:** Tek boyutlu Gauss filtresi (ağırlıklar exp(-d²/(2σ²)), toplamı 1) ve medyan filtresi; kenarlarda yansıtma. Pencere 9, Gauss σ = 2.
- **Madde 4, `dftHarmonics`:** DFT'yi doğrudan tanımından hesaplar, `X[k] = sum x[n] exp(-j 2 pi k n / N)` (sin/cos tablosu bir kere hesaplanır), k = 0..N/2 için harmonik genliklerini döndürür (DC ve N/2 harmoniği |X|/N, diğerleri 2|X|/N; böylece genlik-2 sinüs 2 olarak görünür). Faz tahmini: arg X[k] + 90° (sin için). Doğrulama olarak ters DFT yazıldı; geri elde edilen vektör ile orijinal arasındaki en büyük fark 3.1e-15 çıktı.

## 2. Sonuçlar

### Madde 2: gürültü
En büyük gürültü 0.5981 (sınır 0.6), RMS 0.3376 (kuram 0.3464), gürültülü sinyalin SNR'i 12.44 dB.

### Madde 3: filtreleme (oran 0.3, pencere 9)

| Sinyal | MSE | SNR (dB) | 8 Hz bileşen genliği (gerçek 2.0000) |
|---|---|---|---|
| Gürültülü | 0.1140 | 12.44 | 1.9882 |
| Gauss | 0.0324 | 17.90 | 1.8624 |
| Medyan | 0.0650 | 14.88 | 1.8596 |

Gürültü oranına göre MSE (pencere 9):

| Oran | Gürültülü | Gauss | Medyan |
|---|---|---|---|
| 0.1 | 0.0127 | 0.0135 | 0.0221 |
| 0.3 | 0.1140 | 0.0324 | 0.0650 |
| 0.6 | 0.4558 | 0.0929 | 0.1385 |
| 1.0 | 1.2661 | 0.2337 | 0.3160 |

Pencere boyutuna göre (oran 0.3, Gauss σ = pencere/4):

| Pencere | Gauss MSE | Gauss genlik | Medyan MSE | Medyan genlik |
|---|---|---|---|---|
| 3 | 0.0481 | 1.9707 | 0.0643 | 1.9910 |
| 5 | 0.0325 | 1.9420 | 0.0551 | 1.9543 |
| 9 | 0.0345 | 1.8445 | 0.0650 | 1.8596 |
| 15 | 0.0980 | 1.6140 | 0.1375 | 1.6643 |
| 25 | 0.4408 | 1.0998 | 0.7529 | 0.8794 |
| 41 | 1.4077 | 0.3512 | 3.3655 | 0.5069 |

Ek deney, dürtü gürültüsü (elemanların bir kısmı ±1.5·genliğe ayarlandı, pencere 9):

| Bozuk eleman oranı | Gürültülü MSE | Gauss | Medyan |
|---|---|---|---|
| 0.02 | 0.2633 | 0.0511 | 0.0201 |
| 0.05 | 0.5258 | 0.1163 | 0.0290 |
| 0.10 | 0.8165 | 0.1679 | 0.0352 |

### Madde 4: DFT harmonik vektörü (k. harmonik = k Hz)

| Sinyal | Tepe k | Tepe genlik | Faz (derece) | Taban ortalaması | Taban en büyük | Tepe / taban ortalaması |
|---|---|---|---|---|---|---|
| Temiz | 8 | 2.0000 | 30.00 | yaklaşık 2e-16 | yaklaşık 1e-5 altı | yaklaşık 8.7e15 |
| Gürültülü | 8 | 1.9882 | 29.48 | 0.03833 | 0.07965 | 51.9 |
| Gauss sonrası | 8 | 1.8624 | 29.39 | 0.01324 | 0.08083 | 140.7 |
| Medyan sonrası | 8 | 1.8596 | 29.57 | 0.02487 | 0.08017 | 74.8 |

Gürültülü vektörün en güçlü 5 harmoniği: k = 8 (1.9882), k = 2 (0.0796), k = 98 (0.0787), k = 72 (0.0773), k = 92 (0.0760).

## 3. Yorum

**Gürültü ve filtreleme (madde 3).**

- Düzgün (sınırlı) dağılımlı gürültüde Gauss filtresi medyandan yaklaşık iki kat daha düşük hata verdi (MSE 0.0324 ve 0.0650). Gauss doğrusal bir alçak geçiren filtredir; ortalama alarak rastgele gürültüyü bastırır. Medyan doğrusal olmayan bir filtredir ve bu tür gürültüde ortalamanın sağladığı kazancı sağlamaz.
- Filtreler gürültüyle birlikte sinyalin kendisini de zayıflatır: pencere 9'da 8 Hz bileşeni iki filtrede de yaklaşık %7 azaldı (2.0'dan 1.86'ya). Gauss için kendi hesabımız (kesik çekirdekle) yaklaşık 1.87 veriyor, ölçümle uyumlu.
- Gürültü azken filtre fayda yerine zarar verebilir: oran 0.1'de gürültülü sinyalin MSE'si (0.0127) Gauss (0.0135) ve medyan (0.0221) sonrasından daha düşüktü; filtrenin sinyali bozması, söndürdüğü gürültüden fazlaydı. Gürültü büyüdükçe kazanç artar: oran 1.0'da Gauss 5.4 kat, medyan 4 kat iyileştirdi.
- Pencere boyutu sinyalin periyodundan (32 örnek) çok küçük seçilmelidir. Bu denemede en iyi pencere iki filtre için de yaklaşık 5 çıktı. Pencere büyüdükçe sinyal zayıflar ve bozulur: pencere 41'de Gauss genliği 0.35'e düştü, medyanın MSE'si (3.37) sinyalin gücünden (2.0) bile büyük oldu; yani filtre sonucu, sıfır çıkışından da kötüydü.
- Dürtü gürültüsünde sonuç tersine döner: medyan, Gauss'tan 2.5 ile 4.8 kat düşük hata verdi (örneğin %5 bozuk elemanda 0.0290'a karşı 0.1163). Gauss her bozuk değeri komşularına yayar, medyan ise aykırı değeri sıralamada uca iterek atar. Medyan hatası, dürtü oranı 5 kat artınca (0.02'den 0.10'a) sadece 0.0201'den 0.0352'ye çıktı.

**DFT (madde 4).**

- Temiz sinyalin harmonik vektöründe tek tepe vardır: k = 8, genlik 2.0000, faz 30.00°; diğer bütün harmonikler sayısal hata düzeyindedir (yaklaşık 1e-16). Bunun sebebi 8 tam periyodun vektöre sığması, dolayısıyla spektral sızıntının (leakage) olmamasıdır.
- Gürültülü vektörde tepe aynı yerde ve neredeyse aynı büyüklüktedir (genlik 1.9882, yani %0.6 sapma; faz 29.48°). Gürültü ise bütün harmoniklere dağılır ve bir gürültü tabanı oluşturur: ortalama 0.0383. Beyaz gürültü için beklenen ortalama σ·√(π/N) = 0.0374'tür, ölçümle uyumlu. En güçlü gürültü harmonikleri (k = 2, 98, 72, 92) rastgele konumlardadır ve en büyüğü 0.08 civarıdır, tepenin yaklaşık %4'ü.
- Tepe/taban oranı (51.9, yaklaşık 34 dB) zaman uzayındaki SNR'den (12.44 dB) çok yüksektir. Sinyal enerjisi tek harmonikte toplanırken gürültü enerjisi 128 harmonike yayılır; bu yüzden frekans uzayında sinyali gürültüden ayırmak çok daha kolaydır.
- Filtre sonrası spektrumda Gauss gürültü tabanını 0.0383'ten 0.0132'ye (yaklaşık üçte bire), medyan 0.0249'a indirdi; ancak 8 Hz tepesi de %7 küçüldü (1.8624 ve 1.8596). Tepe/taban oranı Gauss'ta 140.7'ye, medyanda 74.8'e yükseldi. Tabanın en büyük değeri üç durumda da yaklaşık 0.08 kaldı; olası açıklama, en güçlü gürültü harmoniğinin düşük frekansta (k = 2) olması ve Gauss'un alçak geçiren yapısı yüzünden geçmesidir. Bu ayrıca kontrol edilmedi.

## 4. Sınırlılıklar

- Her deney tek bir gürültü gerçekleşmesiyle (sabit tohum) yapıldı; tekrarlanan denemelerin ortalaması alınmadı.
- Gürültü düzgün dağılımlıdır (Gauss dağılımlı gürültü denenmedi); dürtü gürültüsü deneyi ödevin istediği bir madde değildir, yorumu zenginleştirmek için eklenmiştir.
- Tek frekans (8 Hz) ve tek örnekleme (256 Hz, 256 örnek) denendi. Frekans bin sınırlarına denk gelmeseydi (örneğin 8.5 Hz) spektral sızıntı olur ve temiz sinyalin spektrumu da tek tepe olmazdı; bu denenmedi.
- DFT doğrudan tanımdan (O(N²)) hesaplandı, hızlı Fourier dönüşümü (FFT) kullanılmadı.
- Grafikler (cikti_sinyal/*.png) rapora yorumlanmadan eklenmiştir; yorumlar sayısal çıktıya dayanır.
