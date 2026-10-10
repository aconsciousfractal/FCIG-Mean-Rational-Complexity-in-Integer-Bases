# Chiusura della famiglia di radicale 70p

Questa prova usa i profili marcati, la riparazione modulo 70 e il
complemento finito del presente componente. I teoremi sulle famiglie
a tre primi e 30p sono dipendenze distinte del confronto globale.

## 1. Enunciato e modello

**CL5 con certificati CL3.** Siano

\[
B=2^a5^b7^cp^d,\qquad p\ge11\text{ primo},\qquad a,b,c,d,N\ge1.
\]

Allora

\[
\Delta(B,N)<\delta:=\frac2{85}<c_0=\Delta(210,1).
\]

La quantificazione riguarda tutte le valutazioni positive e tutte
le lunghezze, senza un tetto implicito a B o N. La norma è quella
infinito; il denominatore è positivo e unitario sull'intero radicale.
Tutti i residui, compreso zero, hanno peso 1/M. Precisamente,

\[
R_{M,q}(s)=\min_{\substack{d>0,\ (d,q)=1\\a\equiv sd\pmod M}}
\max(|a|,d),\qquad
\nu_q(M)=\frac1M\sum_{s=0}^{M-1}\log R_{M,q}(s)-\frac12\log M.
\]

Quando rad(M)=q si ha μ(M)=νq(M). Sempre
Δ(B,N)=μ(B^N)/log B: il denominatore di questa normalizzazione
non è log(B^N). Si usano logaritmi naturali e intervalli razionali
con arrotondamenti esterni per tutti i confronti numerici decisivi.

## 2. Gli ancoraggi marcati 70 sono nuovi e uniformi

**CL5 + CL3 geometrici.** I due profili ricostruiti danno

\[
70\mid M\quad\Longrightarrow\quad
\nu_{70}(M)\le
\frac{151153169873202211930725852899372149548917}
{1083197534374707740536733696000000000000000}
<\frac7{50}, \tag{A1}
\]

\[
4900\mid M\quad\Longrightarrow\quad
\nu_{70}(M)\le
-\frac{10143530697957057900420872148667939379339}
{1083197534374707740536733696000000000000000}
<0. \tag{A2}
\]

La seconda stima è negativa. Nelle soglie principali basta il cap
semplice 0, che conserva positivi i numeratori dei maggioranti e
semplifica il raccordo monotono. Il cap negativo esatto è comunque
conservato fra i risultati geometrici.

Gli indici di trasferimento sono K=70 e K=4900. Il contratto dei
profili copre ogni M=KL, senza imporre gcd(K,L)=1 o un vincolo
ulteriore sul radicale di L. Perciò gli ancoraggi sono applicabili
anche quando M contiene il quarto primo p. Non si identificano
arbitrariamente il marchio 70 e il radicale completo 70p.

`profiles/PROOF.md` conserva il passaggio finito–generale: tutte
le classi HNF primitive e le molteplicità, copertura delle celle
di riduzione, quattro vertici esatti, copertura angolare con
errore pagato, integrazione e cuspide. La sola tabella di valori
non sostituisce tali argomenti continui.

I dati sono in `profiles/results/q70_k1.json` e `q70_k2.json`.

## 3. Penalità rigorosa del quarto primo

**CL5 con tabelle CL3.** Per p>4 primo, p∤70 e p|M, la nuova
riparazione e il nuovo conteggio periodico dimostrano

\[
0\le\nu_{70p}(M)-\nu_{70}(M)
<\Phi\!\left(\frac{24}{7p}\right)
+\frac{51480}{2107\sqrt M}, \tag{PE70}
\]

con Φ(x)=log(1+x)+x log(1+1/x). I dettagli sono nella prova
`REPAIR70_PROOF.md` e nei suoi verificatori.

La riparazione parte da un minimo marcato u di altezza H. Quando
p divide il suo denominatore, u è primitivo nelle coordinate
reticolari e si completa a una base u,v di determinante M. La
tabella sceglie w=ju+iv con 1≤i≤4, |j+iτ|≤1 e denominatore
unitario modulo 70. Poiché p>4 e p non divide il denominatore
di v, w conserva anche il nuovo vincolo modulo p. Ne segue

\[
R_{M,70p}(s)\le H+\frac{4M}{H}.
\]

Il conteggio periodico specifico di 70 è

\[
2TU_{70}(T)+S_{70}(T)
\le\frac67T^2+\frac{1170}{301}T\qquad(T\in\mathbb Z_{\ge0}).
\]

Da esso si ottiene un maggiorante della massa dei minimi che
richiedono la riparazione. Il conteggio conserva il numeratore
zero e permette ripetizioni; non assume equidistribuzione dei
minimi modulo p. L'integrazione del costo logaritmico conclude
PE70. Le costanti di 42p non sono state trasportate invariatamente.

## 4. Tre regioni infinite e il complemento esatto

**CL5 + confronti CL3.** Φ è crescente nel suo argomento.
La penalità decresce quindi con p e M; i numeratori delle tre
stime seguenti sono positivi. I bordi scalari certificati sono:

| Regione | Bordo usato | Ancoraggio | Upper per Δ, arrotondato esternamente |
|---|---|---:|---:|
| p≥1000, ogni N | B=M=70000 | 7/50 | 0,022896790361759 |
| p≥73, N≥2 | B=5110, M=5110² | 0 | 0,023009054157963 |
| p≥73, B≥3000000, ogni N | B=M=3000000 | 7/50 | 0,023185984069338 |

Tutti sono strettamente inferiori a δ. Il bordo 1000 non deve
essere primo: è un minorante scalare. La prima regione effettiva
parte dal primo 1009. La seconda usa 4900|B^N per N≥2; la terza
usa soltanto M=B^N≥B.

Fuori dai piccoli primi 11,…,71 resta esattamente

\[
73\le p<1000,\quad N=1,\quad
B=2^a5^b7^cp^d<3000000,\quad a,b,c,d\ge1. \tag{F}
\]

I buchi fra primi sono controllati: 71 e 73 sono consecutivi,
così come 997 e 1009. Le classi delle regioni possono sovrapporsi;
la copertura è esaustiva e non richiede monotonia delle vere medie.

## 5. Inventario e consumo del complemento F

**CL3.** `code/inventory70.py` costruisce F moltiplicando tutte
le potenze ammesse e confronta il risultato con una scansione
distinta dei multipli di 70, seguita da fattorizzazione del
cofattore. L'audit costruisce inoltre il dominio separatamente.

| Oggetto | Conteggio |
|---|---:|
| Primi nel complemento | 148 |
| Cofattori 2^a5^b7^c pertinenti | 46 |
| Moduli distinti di F | 3518 |
| Moduli conclusi direttamente da PE70 | 2855 |
| Moduli che richiedono testimoni | 663 |
| Residui rappresentati nei 663 moduli | 135215780 |
| Rappresentanti tramite s↔M−s | 67608553 |
| Basi che usano già A2 a N=1 | 163 |

Quando 4900 divide B, l'ancoraggio A2 è disponibile anche a N=1:
il verificatore lo usa in tutte le 163 righe pertinenti.

Per ogni modulo residuo e ogni 0≤s≤M/2, il produttore scrive una
coppia (a,d). Il consumatore, separato dalla ricerca, verifica
d>0, gcd(d,M)=1, a≡sd mod M e l'altezza. La simmetria pesa 2
ogni record, tranne s=0 e s=M/2. Sono controllati dominio,
cardinalità, massa, ordine e fine esatta del flusso.

Un testimone ammissibile dà un upper del minimo; non occorre
dimostrare l'ottimalità dei testimoni proposti. I logaritmi sono
maggiorati con interi diretti e log M è minorato. Se U è la
somma pesata degli upper logaritmici a scala 10^9 e L il lower
di log M alla stessa scala, il numeratore medio è 2U−ML.

**Tutti i 663 numeratori verificati sono negativi.** Ne segue
μ(M)<0 per questo complemento e quindi Δ(B,1)<δ. Nessun modulo
ha richiesto il secondo stencil di ricerca. Il maggiore degli
upper medi è

\[
-\frac{1289403771357}{17885000000000}
<-0.072094144330835,
\]

a M=35770, p=73. Questo risultato di segno riguarda F nella
rotta dei testimoni; non classifica il segno di tutti gli offset.

Il run ha prodotto e consumato 540868424 byte di testimoni,
con C++17, OpenMP e UBSan, senza diagnostici. I grandi flussi
non sono duplicati nella ZIP: codice, input, risultati esatti,
ricevute e manifest consentono di rigenerarli.

## 6. I sedici piccoli primi e il raccordo per ogni N

**CL5 con CL3 geometrici.** I primi sono

\[
11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71.
\]

Per ognuno, il profilo diretto con marchio completo q=70p dà
νq(M)≤Aq per tutti i multipli M di q, con Aq<δ log q.
Essendo q|B^N e log B≥log q, ne segue Δ(B,N)<δ per tutte le
valutazioni e tutte le lunghezze. Non è necessario un iteratore
o un censimento residuo interno a questi piccoli radicali.

Se Aq è negativo, l'upper normalizzato uniforme si prende pari
a 0; in generale si usa max(0,Aq)/log q con un denominatore
logaritmico minorato rigorosamente. Dividere un cap negativo
per il logaritmo della minima base darebbe invece un confronto
non valido al crescere della base. Il nuovo assemblaggio
controlla esplicitamente questo punto.

`profiles/results/small_prime_closure.json` conserva tutti gli
Aq, le rotte, le divisibilità e il maggiore upper normalizzato.
Il caso q770 passa già al primo livello, con Aq<0,110779378866594.

## 7. Maggiorante comune, affinamento e confronto con il candidato

**CL3/CL5.** Il massimo dei maggioranti delle sedici rotte
piccole, delle tre regioni infinite, delle 2855 rotte analitiche
finite e del bound 0 per i 663 moduli negativi dà

\[
\Delta(B,N)\le
\frac{25221644582229033312307}{1072177116035751725286301}
<0.023523766927132<\frac2{85}.
\]

La rotta analitica che determina questo maggiorante è
B=M=189140, p=193. Questo è il maggiorante della versione base.

Il supplemento `SHARPENED_PE.md` sfrutta separatamente le due
coordinate della riparazione, mantenendo invariati profili,
inventario, soglie e assegnazione dei 663 moduli ai testimoni.
Dimostra

\[
\nu_{70p}(M)-\nu_{70}(M)
<\Phi\!\left(\frac{12}{7p}\right)+\frac9{14p}
+\frac{971503}{44247\sqrt M}. \tag{PE70-sharp}
\]

I tre maggioranti alle soglie diventano, nell'ordine di §4,
0.021193904550578, 0.014633752502645 e 0.018327399494273.
Anche tutte le 2855 rotte analitiche finite migliorano strettamente.
Le 16 rotte dei piccoli primi e i 663 bound negativi rimangono
validi. Il massimo esatto dei nuovi maggioranti è quindi

\[
\boxed{\Delta(B,N)\le\Gamma_{70}
=\frac{503673291285648397571314}{23233333961830267469605761}
<0.021678907216379<\frac2{85}.}
\]

Il margine razionale è

\[
\frac2{85}-\Gamma_{70}
=\frac{12645114755641595659688}{6833333518185372785178165}>0.
\]

La rotta che determina Γ70 è B=M=65030, p=929. Si tratta del
massimo dei maggioranti usati nella prova, senza identificare
un valore o un massimizzatore della vera Δ. Il certificato di
assemblaggio e tutte le rotte affinate si trovano rispettivamente
in `results/closure_certificate.json` e `results/sharp_routes.csv`.

Il candidato 210 viene nuovamente ricostruito dalla definizione.
Il denominatore 1 dà altezza ≤105, quindi la scansione di tutti
i denominatori unitari fino a 105 è esaustiva. Per il prodotto P
delle 210 altezze si verifica con interi P^170>210^18690, che
equivale a c0>2/85. Il relativo codice e certificato sono inclusi.

La copertura precedente conclude il teorema. Nessun passaggio
assume convergenza in N, attainment del supremo globale,
dominanza dei primoriali o monotonìa della vera media nel radicale.
L'uguaglianza C*=c0 resta aperta fuori dai settori chiusi.

## 8. Significato e limiti delle verifiche

CL3 designa qui una verifica finita con aritmetica esatta; CL5
designa gli argomenti matematici scritti, compreso il passaggio
alle famiglie infinite. I controlli dei testimoni geometrici
sono esaustivi, ma produttore e consumatore condividono parte
del kernel di costo. Il controllo geometrico alternativo è
stratificato e viene descritto come tale. I controlli dei
risultati compatti non vengono presentati come replay dei raw.

Questa è una continuazione computazionale con audit AI. Non si
afferma una formalizzazione CL6, una revisione umana specialistica
o una certificazione di novità bibliografica.
