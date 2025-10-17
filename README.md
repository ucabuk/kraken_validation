## 06_KO_validation

Çoklu hizalama ve bitki KO'larının validasyonu için kraken kullanıldı. Öncelikle kraken kısa okumalar üzerinde çalıştırıldı.
Run : 01_kraken.sl. 
Input : R1 and R2.fastq
Output : kraken and kraken.report

ardından Viridiplantae'ye ait okumalar her örnek için kraken sonuçlarından çıkarıldı.
Run : 02_kraken_extraction.sl
Input: kraken_merged, fastp files, kraken_report
output: taxa_specific_reads

Çıkarılan read'ler daha önce Lama için oluşturulmuş non-redundant catalog'a salmon yoluyla hizalandı. Ardından ilgili KO'lar eggnog ve
