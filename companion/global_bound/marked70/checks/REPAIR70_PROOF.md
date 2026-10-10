# Riparazione del denominatore per 70p e per ogni nucleo 2rs

**Data e ambito.** 8 ottobre 2026. Si dimostrano una riparazione uniforme
per i nuclei `q=2rs`, con `r,s` primi dispari distinti, e i lemmi R70,
P70, D70, PE70 necessari alla famiglia `rad(B)=70p`, `p>=11` primo.
La chiusura del complemento finito e gli ancoraggi geometrici marcati
sono dipendenze distinte. Queste note non affermano `C*=c0`, né la
chiusura di tutte le famiglie con quattro o più primi. Non rivendicano
novità, priorità, formalizzazione o revisione specialistica umana.

**Provenienza.** La rotta geometrica e il conteggio delle incidenze
adattano R30/P30/D30/PE30 del dossier `Cstar_QUATTRO_PRIMI_30P_2026-10-08`,
§§2–3, e le successive note R42. La specializzazione modulo 70 viene
ricostruita nel codice nuovo. La prova analitica per tutti i nuclei `2rs`
e il suo conteggio senza tabella sono espliciti qui sotto. Le copie delle
fonti lette e i loro hash sono in `imports/` e `PROVENANCE.json`.

**Classi delle affermazioni.** CL5 indica un argomento matematico interno;
CL3 un controllo finito con aritmetica esatta; CL2 una regressione finita
del modello, che non dimostra per campionamento un enunciato universale.
Una dipendenza dagli ancoraggi della componente dei profili viene indicata
espressamente. I decimali non determinano i segni dei confronti.

## 1. Definizioni e normalizzazione

Per `M>=1`, `q>1` e `s mod M`, si definiscono

\[
R_{M,q}(s)=\min_{\substack{d>0,\ \gcd(d,q)=1\\a\equiv sd\pmod M}}
\max(|a|,d),\qquad
\nu_q(M)=\frac1M\sum_{s=0}^{M-1}\ln R_{M,q}(s)-\frac12\ln M.
\]

Il minimo esiste, perché il denominatore `d=1` è ammissibile. Tutti i
residui, incluso zero, hanno peso `1/M`; il denominatore è positivo e il
logaritmo naturale. Se `rad(M)=q`, allora `nu_q(M)=mu(M)`. Nel problema
delle basi si usa sempre

\[
\Delta(B,N)=\frac{\mu(B^N)}{\ln B},\qquad \delta=\frac2{85}.
\]

Il termine centratore contiene `ln M`; il denominatore della Delta
contiene `ln B`, anche quando `N>1`.

## 2. Lemma R(2rs): una riparazione di raggio quattro senza tabelle

### 2.1 Lemma aritmetico

**CL5.** Siano `q=2rs`, con `r,s` primi dispari distinti, `z` intero e
`0<=tau<=1/2`. Esistono interi `1<=i<=4`, `j` tali che

\[
\gcd(j+iz,q)=1,\qquad |j+i\tau|\le1. \tag{R-arit}
\]

**Dimostrazione.** Se `z` è pari, i tre numeri `z-1,2z-1,3z-1` sono
dispari e coprimi a coppie. I primi due e gli ultimi due hanno gcd 1
per sottrazione; per il primo e il terzo, il gcd divide 2 e un numero
dispari. Ciascuno dei due primi `r,s` può dunque dividere al più uno dei
tre numeri. Uno è coprimo a `q`: si sceglie il corrispondente `i` fra
1,2,3 e `j=-1`.

Se `z` è dispari e unitario modulo `q`, si sceglie `(i,j)=(1,0)`.
Se `z` è dispari ma non unitario, almeno uno dei primi `r,s` divide `z`;
quel primo non divide né `2z-1` né `4z-1`. Questi ultimi due numeri sono
dispari e coprimi fra loro, poiché `4z-1-2(2z-1)=1`. L'altro primo
non può dividerli entrambi. Si sceglie pertanto `i=2` oppure `i=4`,
sempre con `j=-1`.

Per `j=-1` e `1<=i<=4`, entrambi gli estremi di `-1+i tau` su
`[0,1/2]` hanno modulo al più 1; lo stesso vale per la scelta `(1,0)`.
Per `-1/2<=tau<=0`, si applica la costruzione a `(-z,-tau)` e si
inverte il segno di `j`. Il gcd e il modulo dell'affine rimangono
invariati. L'eventuale doppia scelta a `tau=0` non crea problemi.

**CL5 — Limite esatto di questa classe di riparazioni.** Il raggio 4
è necessario se si richiede R-arit per ogni `z,tau`. Per il teorema
cinese dei resti, esiste `z` dispari con `z=0 mod r` e `2z=1 mod s`.
A `tau=1/4`, per ognuno degli indici `i=1,2,3` gli unici interi `j`
con `|j+i/4|<=1` sono `-1,0`. Con `j=0`, il numero `iz` è divisibile
per `r`. Con `j=-1` e `i` dispari, `iz-1` è pari; con `i=2` è
divisibile per `s`. Nessuno è unitario modulo `q`.
Questa necessità riguarda l'ansatz R-arit, non l'ottimalità geometrica
assoluta della costante di riparazione.

### 2.2 Riparazione del vettore reticolare

**CL5 — Enunciato.** Siano `q=2rs` come sopra, `p>4` primo,
`p` non divisore di `q`, e `p|M`. Scelto un minimo marcato
`u=(a,d)` di altezza `H=R_(M,q)(s)`, se `p` non divide `d` allora
`R_(M,qp)(s)=H`. Se `p|d`, allora

\[
R_{M,qp}(s)\le H+\frac{4M}{H}. \tag{R(2rs)}
\]

Non è richiesta l'ipotesi `q|M` per questo lemma; essa verrà richiesta
dagli ancoraggi marcati della specifica applicazione.

**Dimostrazione.** Nel reticolo
`Lambda_s(M)={(a,d):a-sd in MZ}`, si scrive `a=sd+Mk`.
Se `gcd(k,d)>1`, dividere entrambi i coefficienti reticolari per tale
gcd produce un vettore nel medesimo reticolo, con denominatore positivo
ancora coprimo a `q` e altezza minore. Quindi `gcd(k,d)=1`.
La primitività riguarda le coordinate reticolari `(k,d)`; non autorizza
la divisione arbitraria del gcd delle coordinate ambientali `(a,d)`.

Scelti `x,y` con `kx+dy=1`, si pone `v=(b,e)=(sx-My,x)`.
Allora `det(u,v)=M`, e `gcd(d,e)=1`. Sottraendo da `v` un multiplo
intero di `u` si ottiene

\[
v=\tau u+\frac{M}{a^2+d^2}u^\perp,
\qquad u^\perp=(-d,a),\qquad |\tau|\le\frac12.
\]

Si definisce `z=e d^{-1} mod q`, perché `gcd(d,q)=1`, e si applica
R-arit. Per `w=ju+iv`, il denominatore `jd+ie` è unitario modulo `q`.
Se `p|d`, da `p|M` e `a=sd mod M` segue `p|a`; inoltre `p` non
divide `e` perché `gcd(d,e)=1`. Modulo `p`, il nuovo denominatore è
`ie`, diverso da zero dato che `1<=i<=4<p`. Esso è quindi non nullo
e, dopo un eventuale cambio simultaneo di segno, positivo e coprimo a
`qp`. Infine,

\[
\|w\|_\infty
\le |j+i\tau|H+\frac{iMH}{a^2+d^2}
\le H+\frac{4M}{H}.
\]

Si è usato `a^2+d^2>=H^2>0`. Il caso `a=0` è incluso senza alcuna
divisione per `a`. Se invece `p` non divide il denominatore del minimo
scelto, tale minimo è già ammissibile per `qp`; la monotonia della
restrizione dà l'uguaglianza dei due minimi. Non si afferma che il
vettore riparato sia un minimo.

## 3. R70: tabella ricostruita e verificata

R(2rs) applicato a `q=70=2*5*7` dà R70. Per comodità riproducibile,
su `0<=tau<=1/2` l'algoritmo del §2 produce la tabella seguente.

| Classi di z modulo 70 | i | j |
|---|---:|---:|
| 1,3,9,11,13,17,19,23,27,29,31,33,37,39,41,43,47,51,53,57,59,61,67,69 | 1 | 0 |
| 0,2,4,10,12,14,18,20,24,28,30,32,34,38,40,42,44,48,52,54,58,60,62,68 | 1 | -1 |
| 5,6,7,15,16,21,22,26,35,36,45,49,50,55,56,64,65,66 | 2 | -1 |
| 8,46 | 3 | -1 |
| 25,63 | 4 | -1 |

**CL3.** `verify_repair70.py` ricostruisce tutte le 70 righe, verifica
unitarietà e modulo dell'affine agli estremi e controlla anche la
riflessione. Si coprono 140 semintervalli chiusi. A `z=25` oppure `63`
e `tau=1/4`, il certificato enumera l'ostacolo al raggio 3.
R70 dispone inoltre della prova analitica del §2; non dipende soltanto
dalla ricerca della tabella.

## 4. P70: conteggio periodico esatto

Per ogni intero `A>=0`, siano

\[
U(A)=\#\{1\le e\le A:\gcd(e,70)=1\},\qquad
S(A)=\sum_{\substack{1\le e\le A\\\gcd(e,70)=1}}e.
\]

**CL5 con certificato CL3 — Enunciato.** Vale

\[
2AU(A)+S(A)\le \frac67 A^2+\frac{1170}{301}A. \tag{P70}
\]

**Dimostrazione.** La densità delle unità è `rho=24/70=12/35`.
Per `A=70n+j`, `0<=j<70`, le 24 unità nel periodo hanno somma 840 e

\[
U(A)=24n+U(j),\qquad
S(A)=840n^2+70nU(j)+S(j).
\]

Posti `epsilon_j=U(j)-(12/35)j` e
`D_j=2jU(j)+S(j)-(6/7)j^2`, si ottiene l'identità esatta

\[
2AU(A)+S(A)-\frac67 A^2=210\epsilon_j n+D_j.
\]

Le 70 righe del certificato verificano `3epsilon_j<=1170/301` e
`D_j<=(1170/301)j`. Moltiplicare il primo confronto per `70n` e
sommarlo al secondo dimostra P70 per ogni `A`. Il massimo del primo
coefficiente è `132/35`; il secondo raggiunge `1170/301` a `j=43`,
dove `U(43)=16`, `S(43)=376`, `D_43=1170/7`. La costante lineare è
quindi la minima compatibile con questa precisa disuguaglianza
periodica quadratica. A `A=0` entrambi i membri sono zero.

**CL2.** Il controllo diretto di tutti gli interi `0<=A<=7000` è una
regressione distinta. Il passaggio a tutti gli interi usa l'identità
in `n` e i 70 confronti esatti, non l'estrapolazione di un campione.

## 5. D70: densità dei denominatori da riparare

Si sceglie deterministicamente un minimo marcato 70 `(a_s,d_s)` per
ogni residuo, di altezza `H_s`, e si pone `X_s=H_s^2/M`.
Il vincolo di denominatore più forte è sempre il marchio `70p`.

**CL5 — Enunciato.** Per `p>4` primo, `p` non divisore di 70,
`p|M` e ogni reale `t>0`,

\[
F_{\rm bad}(t):=\frac1M\#\{s:p\mid d_s,\ X_s\le t\}
\le\min\left(1,\frac{6t}{7p}
+\frac{1170\sqrt t}{301\sqrt M}\right). \tag{D70}
\]

**Dimostrazione.** Si pongono `T=sqrt(Mt)` e `A=floor(T/p)`.
Un minimo problematico è `(a,d)=(pb,pe)` con `|b|<=A`,
`1<=e<=A` e `gcd(e,70)=1`. Scritti `m=M/p` e `g=gcd(e,m)`,
la congruenza `se=b mod m` ha `pg` soluzioni modulo `M` se `g|b`,
e nessuna altrimenti. Il numero di tutte le incidenze è pertanto al più

\[
p\sum_{\substack{e\le A\\\gcd(e,70)=1}}
g\bigl(2\lfloor A/g\rfloor+1\bigr)
\le p\{2AU(A)+S(A)\}.
\]

Nell'ultimo passaggio si usa `g<=e`. Il termine `+1` conserva il
numeratore zero. Le ripetizioni degli stati fra più testimoni vengono
ammesse, dunque il conteggio è un maggiorante. Per P70 e `A<=T/p`,
esso è al più `6T^2/(7p)+1170T/301`. Dividere per `M` e limitare
con 1 conclude. Se `A=0`, non ci sono denominatori positivi nel
conteggio. Non si assume equidistribuzione dei minimi modulo `p`.

## 6. PE70: costo medio dell'inserimento

Si definisce, per `x>0`,

\[
\Phi(x)=\ln(1+x)+x\ln(1+1/x).
\]

**CL5 — Enunciato.** Sotto le ipotesi di D70,

\[
0\le\nu_{70p}(M)-\nu_{70}(M)
\le\Phi\!\left(\frac{24}{7p}\right)
+\frac{2340\pi}{301\sqrt M}
<\Phi\!\left(\frac{24}{7p}\right)
+\frac{51480}{2107\sqrt M}. \tag{PE70}
\]

**Dimostrazione.** Il vincolo più forte aumenta o conserva ogni
minimo. Se il minimo marcato scelto è già unitario modulo `p`, il
costo è zero. Altrimenti R70 maggiora il costo logaritmico con
`ln(1+4/X_s)`. Poiché la derivata negativa di `ln(1+4/t)` è
`4/[t(t+4)]`, l'integrazione degli indicatori dà

\[
\nu_{70p}(M)-\nu_{70}(M)
\le\int_0^\infty F_{\rm bad}(t)\frac4{t(t+4)}\,dt.
\]

Il centratore `-(1/2)ln M` si cancella esattamente. Per
`A0=6/(7p)` e `B0=1170/(301sqrt(M))`, si applica
`min(1,A0*t+B0*sqrt(t))<=min(1,A0*t)+B0*sqrt(t)`.
Il primo integrale, diviso a `t=1/A0`, vale `Phi(4A0)`;
il secondo vale `2pi*B0`, con il cambio `t=4v^2`.
Entrambi convergono a zero e all'infinito.
Infine,

\[
\frac{22}{7}-\pi
=\int_0^1\frac{x^4(1-x)^4}{1+x^2}\,dx>0
\]

dà il coefficiente razionale esibito. La funzione
`Phi'(x)=ln(1+1/x)>0` mostra che il bound diminuisce separatamente
quando crescono `p` e `M`. Non si usa una monotonia delle vere medie
né un limite asintotico in `N`.

## 7. PE(2rs) senza tabelle periodiche

Questa sezione fornisce un bound uniforme per tutti i nuclei della
forma `q=2rs`, a costo di un coefficiente peggiore di P70.

**CL5.** Per ogni reale `t>=0` e intero dispari positivo `d`, sia

\[
O_d(t)=\#\{m\ge1:m\text{ dispari},\ dm\le t\}
=\left\lfloor\frac{t}{2d}+\frac12\right\rfloor.
\]

Vale `|O_d(t)-t/(2d)|<=1/2`. Per inclusione ed esclusione,

\[
U_q(t)=O_1(t)-O_r(t)-O_s(t)+O_{rs}(t),\quad
\rho=\frac12\left(1-\frac1r\right)\left(1-\frac1s\right),
\]

quindi `|U_q(t)-rho*t|<=2`, per tutti i reali `t>=0`.
Per ogni intero `A>=0`, l'identità della funzione a gradini è

\[
S_q(A)=A U_q(A)-\int_0^A U_q(t)\,dt.
\]

Gli eventuali salti agli estremi non modificano l'integrale. Scrivendo
`U_q(t)=rho*t+E(t)`, `|E(t)|<=2`, si ottiene

\[
2AU_q(A)+S_q(A)
=\frac52\rho A^2+3A E(A)-\int_0^A E(t)\,dt
\le\frac52\rho A^2+8A. \tag{P(2rs)}
\]

La prova delle incidenze del §5 è invariata, usando i coefficienti
`(5/2)rho` e 8. R(2rs) e l'integrazione del §6 danno pertanto, per
ogni primo `p>4`, `p` non divisore di `q`, e `p|M`,

\[
0\le\nu_{qp}(M)-\nu_q(M)
\le\Phi\!\left(\frac{10\rho}{p}\right)
+\frac{16\pi}{\sqrt M}
<\Phi\!\left(\frac5p\right)+\frac{352}{7\sqrt M}. \tag{PE(2rs)}
\]

Si è usato `rho<1/2`. Questo è un enunciato per tutti i nuclei `2rs`,
dimostrato analiticamente. La piccola banca di otto nuclei del codice
è soltanto una regressione: non sostiene da sola il quantificatore
universale. Nessun bound sulla media di partenza `nu_q(M)` segue
automaticamente da PE(2rs). In particolare, non si può iterare
inalterata questa riparazione quando il marchio contiene già quattro
primi distinti.

**Specializzazione ulteriore certificata.** Il calcolo periodico
generico del verificatore dà per `q=66` i coefficienti
`rho=10/33`, `A66=25/33`, `C66=5645/1617`. Dunque lo stesso argomento
prova

\[
\nu_{66p}(M)-\nu_{66}(M)
<\Phi\!\left(\frac{100}{33p}\right)
+\frac{248380}{11319\sqrt M},
\]

sempre per `p>4`, `p` non divisore di 66 e `p|M`.
Non si ricava qui una soglia di esclusione per 66p senza ancoraggi
marcati 66 certificati separatamente.

## 8. Interfaccia con gli ancoraggi e i raccordi della famiglia 70p

Gli ancoraggi necessari alla partizione principale hanno la forma

\[
70\mid M\Longrightarrow\nu_{70}(M)<\frac7{50},\qquad
70^2\mid M\Longrightarrow\nu_{70}(M)<0. \tag{A70}
\]

Questi sono input della componente dei profili: non sono conclusioni della
tabella di riparazione o delle regressioni finite. La loro accettazione
e i confronti razionali dei gate sono registrati separatamente in
`GATES70.md` e `gates70_certificate.json`, prodotti da
`verify_gates70.py` con i due JSON di ancoraggio conservati in `imports/`.

Per `B=2^a5^b7^c p^d`, `a,b,c,d,N>=1`, si ha
`M=B^N`, `rad(M)=70p`, `B>=70p`. Se `N>=2`, allora `70^2|M`,
anche quando le valutazioni di `B` sono diseguali. Più in generale
vale `70^h|M`, dove `h=N*min(a,b,c)`. Per applicare PE70 si usa
prima `nu70`, quindi si arriva a `nu_(70p)=mu`; i due marchi non
vengono identificati prima della riparazione.

Una volta pagati A70 e i confronti del certificatore dei gate, le tre
condizioni sufficienti fissate per la partizione sono:

| Condizione | B minimo | M minimo | Ancoraggio usato |
|---|---:|---:|---:|
| p>=1000, N>=1 | 70000 | 70000 | 7/50 |
| p>=73, N>=2 | 5110 | 26112100 | 0 |
| p>=73, B>=3000000, N>=1 | 3000000 | 3000000 | 7/50 |

Il numeratore di ciascun bound è positivo e decresce in `p,M`;
`ln B` è positivo e cresce in `B`. Il confronto al bordo è pertanto
sufficiente per tutto il dominio della riga. La soglia 1000 è scalare;
non viene considerata un primo. I primi adiacenti sono 997 e 1009.
La seconda soglia è il primo 73, successivo a 71.

Questi raccordi lasciano il complemento `N=1`, `73<=p<1000`,
`B<3000000`, e i primi `11<=p<=71`; la sua chiusura appartiene alle
componenti di profili pieni e di verifica finita. Non è un risultato del
solo PE70.

## 9. Verifiche del modello e limiti del certificato

**CL2/CL3.** Il verificatore principale esamina direttamente la
definizione su 22 coppie `(M,p)` e 5328 istanze `(M,p,s)`. Alcuni
moduli ricorrono con primi diversi: non sono 5328 residui globalmente
distinti. Sono stati riparati 548 minimi problematici, controllando
primitività reticolare, determinante, riduzione di `tau`, congruenza,
unitarietà, segno e bound sull'altezza. I 32 punti di discontinuità
delle densità problematiche sono confrontati esattamente con D70.
La tabella, i periodi e questi controlli producono 47435 verifiche
aritmetiche nel produttore corrente.

**CL3 — Separazione dei modelli.** Per `M=770,s=70`, il minimo
marcato 70 è `(a,d)=(0,11)`, altezza 11, mentre il minimo pieno ha
altezza 70. Dividere il vettore marcato per 11 produce `(0,1)`, che
non appartiene al reticolo. Questo esempio controlla due errori
possibili: identificare il marchio rilassato con il modello pieno,
oppure dividere un gcd ambientale senza controllare la congruenza.

`verify_repair70.py` è autonomo e usa soltanto la libreria standard
Python. Gli helper aritmetici adattano il codice precedente R42;
il suo replay non costituisce una seconda implementazione indipendente.
Nessuna di
queste verifiche finite viene presentata come sostituto della prova
analitica, di una formalizzazione o della review specialistica umana.
