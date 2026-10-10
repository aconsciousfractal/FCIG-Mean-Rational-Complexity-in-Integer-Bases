# Tre raccordi esatti per la famiglia 70p

**Stato: PASS, 8 ottobre 2026.** `verify_gates70.py` ricertifica i tre
raccordi fissati per la partizione della famiglia. La prova dell'inserimento
è `REPAIR70_PROOF.md`, §§2–6. Gli ancoraggi marcati provengono dalla
ricostruzione nuova della componente dei profili, con produttore, consumatore
esaustivo e audit secondario dichiarati concordi; questo componente verifica
il loro uso aritmetico e conserva le copie con hash. Non riesegue le
geometrie dei profili.

## 1. Ancoraggi importati nella loro forma uniforme

La fonte `imports/q70_k1.json`, per `70|M`, fornisce

\[
\nu_{70}(M)\le
\frac{151153169873202211930725852899372149548917}
{1083197534374707740536733696000000000000000}
<\frac7{50}.
\]

La fonte `imports/q70_k2.json`, per `4900|M`, fornisce

\[
\nu_{70}(M)\le
-\frac{10143530697957057900420872148667939379339}
{1083197534374707740536733696000000000000000}
<0.
\]

Il certificatore controlla i campi `q=70`, `k`, `transfer_index=70^k`
e i due confronti fra frazioni; ne registra byte e SHA-256. Il solo
controllo di questi campi non prova la geometria o la formula di
trasferimento, che restano nella distinta componente dei profili. Non vengono
usati gli ancoraggi storici orientativi `0.22` e `0.09`.

## 2. Bound e quantificatori

Si considerino

\[
B=2^a5^b7^cp^d,\qquad
a,b,c,d,N\ge1,\qquad p\ge11\text{ primo},\qquad M=B^N.
\]

Per PE70, se un ancoraggio applicabile vale `nu70(M)<A`,

\[
\Delta(B,N)
<\frac{A+\Phi(24/(7p))+51480/(2107\sqrt M)}{\ln B},
\qquad
\Phi(x)=\ln(1+x)+x\ln(1+1/x).
\]

Il termine di partenza è una media marcata. Solo dopo l'inserimento
si identifica `nu_(70p)(M)=mu(M)`, perché `rad(M)=70p`.
Il denominatore è `ln B` per tutti gli `N`.

Con gli ancoraggi semplici `A=7/50` e `A=0`, tutti i numeratori sono
positivi. Il numeratore decresce separatamente in `p,M`, e il
denominatore positivo cresce in `B`. È pertanto sufficiente verificare
un bordo inferiore di ciascuna variabile nelle tre regioni seguenti.

| Regione sufficiente | p al bordo | B minimo | M minimo | A |
|---|---:|---:|---:|---:|
| p>=1000, N>=1 | 1000 | 70000 | 70000 | 7/50 |
| p>=73, N>=2 | 73 | 5110 | 26112100 | 0 |
| p>=73, B>=3000000, N>=1 | 73 | 3000000 | 3000000 | 7/50 |

Per la seconda riga, `4900|B^N` segue da `N>=2`; non si impone che le
valutazioni siano uguali. La soglia 1000 è uno scalare del confronto,
non un numero supposto primo o una base concreta della famiglia.

## 3. Certificati numerici esterni

Siano `q0=floor(sqrt(M_min))` e `L0` un minorante razionale positivo
di `ln(B_min)`. Ogni riga è verificata calcolando esattamente

\[
U=\frac{A+\Phi_{\rm upper}(24/(7p_{\min}))
+51480/(2107q_0)}{L_0}.
\]

Il codice verifica `q0^2<=M_min<(q0+1)^2`, e `U<2/85` mediante un
confronto fra interi. Gli upper e i margini esatti completi sono nel
JSON; qui si riportano solo arrotondamenti nella direzione sicura.

| Regione | q0 | Upper decimale esterno per Delta | Minorante decimale del margine da 2/85 |
|---|---:|---:|---:|
| p>=1000, N>=1 | 264 | 0.022896790361759 | 0.000632621402947 |
| p>=73, N>=2 | 5110 | 0.023009054157963 | 0.000520357606743 |
| p>=73, B>=3000000, N>=1 | 1732 | 0.023185984069338 | 0.000343427695368 |

Gli upper della tabella sono maggioranti, non valori effettivamente
assunti da Delta. Il loro massimo riguarda solo queste tre rotte,
non la famiglia intera né il suo complemento finito.

Gli intervalli logaritmici sono razionali: dopo riduzione binaria
`x=2^k y`, `1<=y<2`, si usa

\[
\ln y=2\sum_{j=0}^{T-1}\frac{z^{2j+1}}{2j+1}+E_T,
\quad z=\frac{y-1}{y+1},\quad
0\le E_T\le\frac{2z^{2T+1}}{(2T+1)(1-z^2)},\quad T=80.
\]

La medesima serie con `z=1/3` racchiude `ln2`; per `k<0` si
scambiano correttamente gli estremi del contributo `k ln2`.
Gli intervalli di Phi si sommano con coefficienti positivi.
Il floating point non decide nessun gate.

## 4. Partizione residua e controllo negativo mirato

I primi adiacenti alle soglie sono esattamente `71,73` e `997,1009`,
come verifica il codice per divisione intera. Rimangono da trattare
con le altri componenti:

1. I 16 primi `11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71`,
   per tutti gli esponenti e tutte le lunghezze.
2. `N=1`, i 148 primi `73<=p<1000`, e le basi della famiglia con
   `B<3000000`.

Questa è una partizione finita dopo aver separato i piccoli primi;
nessun gate pretende di avere già esaminato le basi residue.

Il codice certifica inoltre che il medesimo bound razionale sufficiente
con ancoraggio semplificato **0**, al bordo `p=71`, `B=4970`, `N=2`,
è maggiore di `2/85`. Si tratta soltanto di una limitazione di quel
bound semplificato. Non è un controesempio per Delta e non esclude
che l'ancoraggio negativo più preciso o un profilo pieno chiuda il caso.

## 5. Riproduzione e portabilità

Nella directory contenente i file:

```bash
python verify_repair70.py
python verify_gates70.py
```

Entrambi usano soltanto la libreria standard Python e scrivono i JSON
nella propria directory. Il secondo importa il primo come modulo e
richiede la sottodirectory `imports/` adiacente. Il manifest conserva
gli hash dei sorgenti e dei risultati finali. Il controllo rapido non afferma revisione specialistica indipendente.
