# Coda uniforme a cinque primi: rad(B)=770 ell, ell>=300

**Stato: PASS, 8 ottobre 2026.** Con una nuova tabella modulo 770,
un nuovo conteggio periodico e il profilo marcato 770 del componente
profili70, si dimostra il risultato seguente.

## 1. Enunciato e dipendenze

**CL5 con certificati CL3 e ancoraggio distinto.** Per ogni primo
`ell>=300` e per tutti gli interi positivi `a,b,c,d,e,N`,

\[
B=2^a5^b7^c11^d\ell^e
\quad\Longrightarrow\quad
\Delta(B,N)<\frac2{85}. \tag{F770}
\]

300 è una soglia scalare, con primo effettivo successivo 307.
Le valutazioni possono essere arbitrarie e diseguali. Il teorema
richiede il contratto marcato

\[
770\mid M\quad\Longrightarrow\quad\nu_{770}(M)<\frac9{80}
\]

per **tutti** i multipli `M` di 770, compresi quelli con il primo
aggiuntivo `ell`. Questo è un ancoraggio della media rilassata
`nu770`, non una mera riformulazione del teorema sulle basi con
radicale esattamente 770.

Si usano le definizioni di `REPAIR70_PROOF.md`, §1, in particolare
`Delta(B,N)=mu(B^N)/ln B`. La normalizzazione non cambia per `N>1`.
Questa nota non chiude i radicali con cinque primi in generale,
né tutti i primi `ell` della forma `770ell`; non determina `C*`.
Il confronto distinto `2/85<c0` appartiene al certificato del
candidato della componente completo.

## 2. R770: conservare il marchio dei quattro primi precedenti

**CL5 con tabella CL3.** Siano `ell>6` primo, `ell` non divisore
di 770 e `ell|M`. Per un minimo marcato 770 di altezza `H`, se il
denominatore è divisibile per `ell`, si ha

\[
R_{M,770\ell}(s)\le H+\frac{6M}{H}. \tag{R770}
\]

Se il denominatore del minimo scelto è già coprimo a `ell`, i due
minimi coincidono. L'ipotesi `770|M` non è necessaria a R770.

**Dimostrazione.** L'argomento delle coordinate reticolari del §2.2
della prova70 vale per ogni marchio: il minimo `u=(a,d)` si completa
a una base `u,v=(b,e)` con determinante `M`, `gcd(d,e)=1` e

\[
v=\tau u+\frac{M}{a^2+d^2}u^\perp,\qquad |\tau|\le\frac12.
\]

Per `z=e d^{-1} mod770`, serve una tabella **modulo 770**. La
tabella70 non conserva da sola il primo 11. La nuova tabella copre
ogni `z` e ogni `tau` con una coppia `1<=i<=6`, `j` tale che
`gcd(j+iz,770)=1` e `|j+i tau|<=1`.

Per `0<=tau<=1/2`, nelle classi ordinarie si sceglie `(1,0)` quando
`z` è un'unità, altrimenti il primo `i` fra 1,2,3,4 con `iz-1`
unitario e `j=-1`. Le sei classi speciali sono coperte da:

| Classi z modulo 770 | Intervallo chiuso di tau | (i,j) |
|---|---|---|
| 25,259,319,553,583,765 | [0,1/3] | (6,-1) |
| 25,553,583,765 | [1/3,1/2] | (3,-2) |
| 259,319 | [1/5,1/2] | (5,-2) |

L'ultima riga si sovrappone alla prima fra `1/5` e `1/3`, senza
lasciare buchi. Per `tau<0`, si riflettono `z,tau` e si cambia il
segno di `j`. Lo zero e tutte le giunzioni sono inclusi.

**CL3.** `verify_repair770.py` ricostruisce le coppie da intervalli
razionali, poi ricontrolla ogni gcd, gli estremi affini, la copertura
e la riflessione. Le 770 classi producono 776 intervalli positivi,
1552 con riflessione. Il certificato conserva tutte le righe.

**CL5.** Per `w=ju+iv`, il denominatore è unitario modulo 770.
Poiché `ell|d` e `gcd(d,e)=1`, `ell` non divide `e`; pertanto
`jd+ie=ie mod ell` è non nullo, perché `1<=i<=6<ell`.
Un eventuale cambio simultaneo di segno rende positivo il denominatore.
Come nella prova70,

\[
\|w\|_\infty
\le |j+i\tau|H+\frac{iMH}{a^2+d^2}
\le H+\frac{6M}{H}.
\]

Il caso `a=0` è incluso; non viene fatta alcuna divisione per `a`
o per un gcd ambientale non giustificato.

**CL3 — Limite dell'ansatz.** A `z=25`, `tau=1/10`, con
`1<=i<=5` gli unici `j` possibili sono `-1,0`. Con `j=0`, `iz`
è divisibile per 5. Con `j=-1`, gli indici dispari producono un
numero pari, `i=2` produce 49 e `i=4` produce 99. Nessuno è unitario
modulo 770. Il codice verifica questi dieci casi. Il raggio 6 è
quindi necessario per questa classe di tabelle; non si afferma
l'ottimalità geometrica assoluta del costo.

## 3. P770, D770 e PE770

**CL5 con 770 confronti CL3.** Le 240 unità modulo 770 hanno densità
`rho=24/77` e somma per periodo 92400. Per `U,S` definiti come
nel conteggio P70, ma con marchio 770, vale per ogni intero `A>=0`

\[
2AU(A)+S(A)\le\frac{60}{77}A^2+\frac{453}{77}A. \tag{P770}
\]

Infatti, per `A=770n+j`, `0<=j<770`,

\[
U(A)=240n+U(j),\qquad
S(A)=92400n^2+770nU(j)+S(j).
\]

Posti `epsilon_j=U(j)-(24/77)j` e
`D_j=2jU(j)+S(j)-(60/77)j^2`, il residuo dal termine quadratico
è `2310epsilon_j n+D_j`. Il certificato verifica
`3epsilon_j<=453/77` e `D_j<=(453/77)j`. La prima costante
è raggiunta a `j=183`. Moltiplicare il primo confronto per `770n`
e sommare il secondo dimostra P770 per ogni periodo.

Il medesimo conteggio delle incidenze di D70 conserva il numeratore
zero e ammette eventuali duplicati. Con `X=H^2/M`, dà

\[
\Pr(\ell\mid d,\ X\le t)
\le\min\left(1,\frac{60t}{77\ell}
+\frac{453\sqrt t}{77\sqrt M}\right). \tag{D770}
\]

Si usano soltanto la congruenza ridotta modulo `M/ell` e
`gcd(e,M/ell)<=e`; non si assume equidistribuzione.

R770 limita il costo logaritmico con `ln(1+6/X)`. Il suo nucleo
di integrazione è `6/[t(t+6)]`. Integrare il termine lineare
troncato di D770 dà `Phi(360/(77ell))`; integrare il termine
radicale dà `pi*sqrt(6)*453/(77sqrt(M))`. Quindi

\[
\begin{split}
0\le\nu_{770\ell}(M)-\nu_{770}(M)
&\le\Phi\!\left(\frac{360}{77\ell}\right)
+\frac{453\pi\sqrt6}{77\sqrt M}\\
&<\Phi\!\left(\frac{360}{77\ell}\right)
+\frac{453}{10\sqrt M}. \tag{PE770}
\end{split}
\]

Il centratore `-(1/2)ln M` si cancella fra i due marchi sullo stesso
modulo. Il maggiorante razionale usa `pi<22/7`, dimostrato
nell'integrale positivo della prova70, e `sqrt6<49/20`, dato che
`2401>2400`. Esattamente,

\[
\frac{22}{7}\frac{49}{20}=\frac{77}{10},\qquad
\frac{453}{77}\frac{77}{10}=\frac{453}{10}.
\]

## 4. Ancoraggio marcato e coda scalare 300

La fonte congelata `imports/q770.json` è il nuovo profilo marcato
770, orizzonte 1, della componente dei profili. Il suo `total_upper` esatto è

\[
A_{770}=
\frac{119995950047855738283362404860350854113141}
{1083197534374707740536733696000000000000000}
<\frac9{80}.
\]

Il contratto geometrico trasferisce questo bound a `nu770(M)` per
ogni `770|M`, anche se `M` possiede altri divisori primi. Il nuovo
codice verifica l'identità del marchio, l'indice di trasferimento,
la frazione e la sua semplificazione a `9/80`, conservando byte e hash.
La validità geometrica uniforme è una dipendenza della componente dei profili;
non viene dedotta dal solo JSON, né dalla sola chiusura a quattro
primi con radicale 770.

Per la famiglia F770, `B>=770ell`, `M=B^N>=B` e `770|M`.
Il marchio finale `770ell` coincide con `rad(M)`. Perciò

\[
\Delta(B,N)
<\frac{\frac9{80}+\Phi(360/(77\ell))+
\frac{453}{10\sqrt M}}{\ln B}.
\]

Il numeratore è positivo e decresce separatamente in `ell,M`;
`ln B` è positivo e cresce in `B`. Basta il bordo scalare
`ell=300`, `B=M=231000`. Con

\[
480^2=230400\le231000<231361=481^2,
\]

il verificatore costruisce un upper interamente razionale del rapporto,
usando `sqrt(M)>=480` e gli intervalli logaritmici descritti in
`GATES70.md`, §3. Ottiene

\[
\Delta(B,N)<0.023273724251698<\frac2{85},
\]

con margine razionale maggiore di `0.000255687513008` da `2/85`.
Entrambi i decimali sono arrotondati verso l'esterno; l'upper e il
margine completi sono nel JSON. L'upper non è un valore effettivo
di Delta. La verifica dei primi mostra che 293 e 307 sono adiacenti;
il primo realmente coperto oltre lo scalare 300 è 307.

La monotonia del bound copre tutti gli esponenti positivi e tutte
le lunghezze. Restano fuori da questa nota gli altri primi `ell`
non divisori di 770 sotto 300, incluso 3, e le altre famiglie a
cinque primi. Il lemma di inserimento con raggio 6 richiede inoltre
`ell>6`; non viene applicato al primo 3.

## 5. Regressioni finite e riproduzione

Il verificatore controlla direttamente il modello su tre coppie
`(M,ell)=(10010,13),(13090,17),(169,13)`: 23269 istanze `(M,ell,s)`,
1160 minimi riparati e 58 punti delle densità. Il modulo 169 mostra
anche l'applicazione di R770 senza `770|M`; non è usato come base
della famiglia F770. Il controllo periodico diretto fino ad `A=7700`
è una regressione distinta dal passaggio rigoroso a tutti i periodi.

Il codice adatta il produttore precedente, e importa gli helper locali di `verify_repair70.py`
e `verify_gates70.py`. Non è un audit indipendente della propria
implementazione. Il replay non afferma revisione specialistica indipendente.

```bash
python verify_repair770.py
```

Lo script usa soltanto la libreria standard Python, legge
`imports/q770.json` adiacente e scrive `repair770_certificate.json`
nella propria directory. La prova, la tabella finita, i test del modello
e la dipendenza dagli ancoraggi rimangono distinti.
