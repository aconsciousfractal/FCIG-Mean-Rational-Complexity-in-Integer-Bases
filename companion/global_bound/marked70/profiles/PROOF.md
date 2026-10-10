# Profili marcati 70 e piccoli radicali 70p

## 1. Dominio e contratto

Per un marchio squarefree q definiamo

\[
R_{M,q}(s)=\min_{d>0,\,(d,q)=1,\ a\equiv sd\pmod M}\max(|a|,d),\qquad
\nu_q(M)=\frac1M\sum_{s=0}^{M-1}\log R_{M,q}(s)-\frac12\log M.
\]

Si includono tutti i residui, compreso zero. I logaritmi sono naturali. Se
q=rad(M), allora νq(M)=μ(M), e Δ(B,N)=μ(B^N)/log B. Il confronto scientifico
di questo componente è δ=2/85.

**CL5 con verifiche finite CL3.** I profili completi coprono esattamente
i sedici primi

`11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71`.

Per ciascuno di essi e per tutti a,b,c,d,N≥1,

\[
\Delta(2^a5^b7^cp^d,N)\le\Gamma<2/85.
\]

Il razionale Γ, il margine e tutte le rotte sono ricostruiti da
`code/assemble_small.py` e registrati in `results/small_prime_closure.json`.
Γ è un maggiorante delle rotte, non un massimo osservato della media.
Le valutazioni dei quattro primi non devono essere uguali.

**CL5 con CL3.** Le due nuove ricostruzioni marcate danno anche

\[
70\mid M\Longrightarrow\nu_{70}(M)\le A_{70,1}<7/50,
\qquad
4900\mid M\Longrightarrow\nu_{70}(M)\le A_{70,2}<0,
\]

con i cap razionali esatti in `results/q70_k1.json` e `results/q70_k2.json`.
Qui M può contenere ulteriori primi: il denominatore è unitario modulo
il marchio 70, non necessariamente modulo rad(M). La riparazione verso
70p appartiene al componente PE70 separato.

Il dominio p≥11 esclude p=3, che darebbe il radicale 210, e i primi
già presenti in 70. Nessun cap storico 30p, 42p o T3 è assunto. Sono
riutilizzati il metodo analitico e il codice aritmetico, con provenienza
esplicita in `imports/provenance.json` e `imports/initial_42_to_70.diff`.

## 2. Trasferimento con il marchio conservato

Per una base per righe C di determinante positivo e una colonna m
primitiva modulo ogni primo di q, poniamo

\[
F_q(C,m)=\min_{(nm,q)=1}\log\frac{\|nC\|_\infty}{\sqrt{\det C}}.
\]

Si usa K=q^k, con k=1 per q=70p e k=1,2 per q=70. Enumeriamo tutte
le HNF H=(a,0;b,d), ad=K, 0≤b<a, non nulle modulo ciascun primo del
marchio. A ogni primo hanno rango uno. Ogni etichetta kernel in
∏ℓ|q P¹(Fℓ) compare K/q volte. Le etichette sono 144 per q=70 e
144(p+1) per q=70p. Pertanto il numero totale dei rami è (K/q) volte
il numero delle etichette; al secondo orizzonte marcato 70 sono 10.080.

Un marchio genitore ammette le etichette diverse dalla propria retta
in ogni fattore primo, ossia esattamente K rami. I marchi immagine di
uno stesso ramo sono multipli unitari e hanno gli stessi testimoni.
La somma sulle etichette è il prodotto tensoriale ∏(Jℓ+1−Iℓ+1), senza
eliminare valori negativi. La molteplicità K/q e la normalizzazione K
sono verificate per tutti i marchi; l'audit Python enumera di nuovo
le rette locali, le immagini e i kernel con un algoritmo distinto.

Sul reticolo fisico C=(L,0;s,1), m=(0,1)ᵀ, l'ammissione impone d=1.
I figli sono quindi precisamente s+bL modulo KL, 0≤b<K. Ne segue

\[
\nu_q(KL)=\frac1L\sum_{s=0}^{L-1}
(T_{q,K}F_q)(C_{L,s},m_0).
\]

Non occorre gcd(L,q)=1. La dimostrazione completa del trasferimento e
dell'invarianza al cambio simultaneo di base e marchio è conservata
in `imports/DIRECT_METHOD_SOURCE.tex`, sezioni “Trasferimento marcato
senza perdita di stati” e “Dalle celle a ogni modulo iniziale”. È un
argomento per marchi squarefree, applicabile ai fattori 2,5,7,p.

## 3. Copertura geometrica continua e consumo dei testimoni

Una base ridotta di Gauss, a meno di scala e rotazione, ha forma

\[
Y(\rho,t)=\begin{pmatrix}\rho&0\\t\rho&1\end{pmatrix},\quad
0<\rho\le2/\sqrt3,\quad |t|\le1/2.
\]

La griglia razionale comprende 128 intervalli radiali e 21 intervalli
di taglio: 2.688 celle, denominatore di taglio T₀=40. `assemble.grid`
verifica tutte le etichette, la continuità radiale e di taglio e
l'estremo superiore di Gauss. La parte sotto il primo raggio è
trattata dalla cuspide, non scartata.

Il produttore propone coefficienti interi per ciascuna coppia
cella/ramo. Il consumatore è compilato con `BUILD_CONSUMER`, senza
il motore di ricerca, e legge tutti i record: lunghezza esatta,
unitarietà del marchio, non annullamento e intervalli aritmetici
sono verificati. Rigenera poi tutti i cap e le somme tensoriali.

Per i coefficienti (u,v) nella base genitrice, al nodo tanθ=n/A,
ρ=x/D₀ e t=j/T₀, calcola

\[
U=Ax(T_0u+jv)+T_0D_0nv,\quad
V=-nx(T_0u+jv)+T_0D_0Av,\quad R=\max(|U|,|V|),
\]

\[
\log Z=\log R-\tfrac12\{\log T_0^2+\log D_0+
\log(A^2+n^2)+\log K+\log x\}.
\]

Il logaritmo di R è maggiorato e quelli sottratti sono minorati;
le divisioni superiori sono corrette anche per interi negativi.
Tutte le decisioni usano interi o razionali con controllo dei range.

A θ fissato, Z è separatamente convessa in t e in z=(logρ)/2:
ogni componente normalizzata ha forma A(t)e^z+Be^−z e il suo valore
assoluto è convesso. Il massimo della norma sulla cella è dunque
raggiunto a un vertice, prima del logaritmo crescente. Il centro
serve solo a proporre i coefficienti: non occorre che la sua
approssimazione sia esatta o interna, perché i quattro vertici sono
certificati nuovamente.

Per ogni vettore f(θ)=log‖vRθ‖∞ soddisfa f″≥−2 nei tratti regolari,
con salti della derivata verso l'alto. La proprietà passa al massimo
sui vertici e sui marchi e alle medie. I nodi θn=arctan(n/A),
−A≤n≤A, coprono quindi tutte le orientazioni pagando 1/(4A²).
Si usa A=4 per tutti i sedici marchi completi e A=16 per i due
profili marcati 70. La norma resta la norma infinito.

## 4. Dal profilo a tutti i moduli

Per i reticoli fisici di determinante L, se λs è il vettore euclideo
nonzero più corto, ρs=λs²/L soddisfa

\[
\Pr(\rho_s\le x)\le\min(1,\alpha x),\qquad\alpha=29/14.
\]

La dimostrazione per incidenze, riportata nel lemma SV della fonte
diretta allegata, non assume indipendenza o equidistribuzione. Per
Y=√(Lx)<L il vettore orizzontale è escluso e la congruenza per ogni
denominatore d ha gcd(d,L) soluzioni quando quel divisore divide il
numeratore. Il conteggio nel disco è al più
2∑d√(Y²−d²)+n(n+1)/2, n=⌊Y⌋. Concavità e trapezi lo maggiorano
con (π+1)Y²/2<29Y²/14. I casi Y<1 e Y≥L sono separati e le
eventuali ripetizioni sono ammesse nella maggiorazione.

Il CRT dà un testimone con coefficienti in {0,…,q−1}. Nei
coefficienti genitori i moduli sono al più 2(q−1)K e (q−1)K.
Per ogni orientazione vale quindi la cuspide

\[
T_{q,K}F_q\le\log(q-1)+\tfrac12\log K-
\tfrac12\log\rho+\tfrac52\rho.
\]

Sia b₀≥⋯≥bm−1 l'inviluppo monotono dei massimi per fascia radiale,
rj gli estremi superiori e ρ₀ il primo estremo inferiore. Si calcolano

\[
\mathcal B=b_{m-1}+\sum_{j<m-1}(b_j-b_{j+1})\min(1,\alpha r_j),
\]
\[
\mathcal E=\alpha\rho_0\left(\max\{c-b_0,0\}+\tfrac12\right),\quad
c=\log(q-1)+\tfrac12\log K-\tfrac12\log\rho_0+\tfrac52\rho_0.
\]

L'integrazione della CDF dà E[log(ρ₀/ρ)₊]≤αρ₀. Pertanto, per
**ogni** L≥1, νq(KL)≤B+E. È questo passaggio analitico continuo,
insieme ai cap uniformi delle celle, a estendere il risultato a tutti
i moduli; i banchi di piccoli moduli sono soltanto controlli ulteriori.

## 5. Chiusura esaustiva dei sedici radicali

Per ogni primo dichiarato il cap diretto Aq,1 è minore di δlogq.
Se B ha radicale q=70p, q divide B^N per ogni N≥1, mentre logB≥logq.
Usando la parte positiva del cap quando necessario si ottiene

\[
\Delta(B,N)\le\frac{\max(0,A_{q,1})}{\log q}<\delta.
\]

Il verificatore controlla tutti i primi da 11 a 71, senza inferire
monotonia rispetto a p. Il massimo di questi sedici maggioranti è Γ.
Non servono uguaglianza delle valutazioni, raffinamenti, profili
iterati a marchio completo o un complemento di basi finite interno.

I due profili del nucleo restano distinti: 70|M è sufficiente per
7/50, mentre il cap 0 richiede 70²|M. In particolare N≥2 implica
questa seconda divisibilità per ogni B di radicale 70p, anche con
valutazioni diseguali. La famiglia p≥73 è demandata a PE70 e al
complemento finito del pacchetto principale.

## 6. Riproducibilità e limiti

Produttore e consumatore devono rigenerare output identici; il
manifest di ciascun profilo registra lunghezza e SHA256 dei testimoni.
Il replay completo rigenera i flussi in una directory temporanea,
li consuma e confronta anche hash, algebra, integrali e audit.
Le tabelle dei logaritmi sono controllate con una serie razionale
indipendente. L'audit Python ricostruisce tutti i rami e tutti i
valori tensoriali costanti; i test dei tensori firmati e della
geometria alternativa sono dichiaratamente stratificati. Per ogni
profilo si ricontrollano 96 testimoni a tutti i nodi e vertici,
oltre a punti razionali interni. I banchi definitori enumerano tutti
gli stati dei piccoli moduli dichiarati nei JSON indipendenti.

I controlli negativi respingono flussi troncati, record extra,
coefficienti nonunitari, un parametro composto e un buco radiale.
Un vettore esplicito dimostra che la griglia angolare finita da sola
non certifica il massimo continuo: il costo angolare è necessario.

Il kernel aritmetico di producer e consumer è condiviso. Non si
afferma una seconda implementazione geometrica esaustiva né una prova
formale CL6 o una revisione specialistica indipendente. I decimali
sono visualizzazioni di razionali; il tempo di esecuzione è l'unico
uso diagnostico dei floating point nel programma geometrico.

La verifica rapida reintegra tutte le celle memorizzate e controlla
il raccordo e i logaritmi, ma non ricostruisce i testimoni geometrici
né riesegue i test negativi. Questo componente non identifica C*.
