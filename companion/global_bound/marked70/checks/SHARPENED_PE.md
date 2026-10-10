# Miglioramento geometrico del costo di inserimento

**Stato: PASS, 8 ottobre 2026.** Il miglioramento sfrutta separatamente
le due coordinate del vettore minimo. Usa le stesse tabelle, gli stessi
conteggi di densità e gli stessi profili già certificati. Non cambia
cutoff o inventari. Può essere usato come secondo maggiorante sulle
medesime rotte della componente completo.

La prova puntuale, l’integrazione e le costanti sono distinte dai controlli aritmetici finiti. Non si afferma una revisione specialistica indipendente.

## 1. Ipotesi precise

Si suppone una tabella di riparazione per un marchio `q` e un raggio
intero `r>=1`, che fornisca, per ogni `z mod q` e `|tau|<=1/2`,
interi `1<=i<=r`, `j` con

\[
\gcd(j+iz,q)=1,\qquad |j+i\tau|\le1.
\]

Il nuovo primo soddisfa `p>r`, `p` non divisore di `q`, `p|M`.
Si usa la struttura del vettore riparato, non soltanto un bound
astratto sulla sua altezza. Per un minimo problematico `u=(a,d)`,
`H=max(|a|,d)`, il vettore prima del cambio simultaneo di segno è

\[
w=\alpha u+\frac{iM}{a^2+d^2}u^\perp,
\qquad \alpha=j+i\tau,\quad |\alpha|\le1.
\]

La tabella garantisce che il denominatore riparato è non nullo e
coprimo a `qp`, come nelle prove R70 e R770. Si assume inoltre il
bound già dimostrato per la densità dei minimi problematici

\[
F_{\rm bad}(X)\le\min(1,A X+B\sqrt X),\qquad A>0,\quad B\ge0.
\]

Nelle applicazioni, `A=A_q/p` e `B=C_q/sqrt(M)`, con i coefficienti
del conteggio periodico. Per il raccordo dell'integrale lineare si
richiede `y=rA<=4/3`. Questa ipotesi è automatica nei domini di
PE70, PE770 e del lemma uniforme per nuclei `2uv`.

## 2. Un maggiorante puntuale più forte

Si pongono

\[
X=H^2/M,\quad v=r/X,\quad
\eta=\min(|a|,d)/H\in[0,1].
\]

Prendendo i moduli delle due coordinate separatamente e usando
`i<=r`, si ha

\[
\frac{\|w\|_\infty}{H}
\le\max\left(1+\frac{v\eta}{1+\eta^2},
\eta+\frac{v}{1+\eta^2}\right). \tag{C}
\]

**CL5.** Il membro destro è al più

\[
G(v)=
\begin{cases}
1+v/2,&0<v\le4/3,\\
\sqrt{1+v^2},&v\ge4/3.
\end{cases} \tag{G}
\]

**Dimostrazione.** La prima coordinata di C è al più `1+v/2`
perché `2eta<=1+eta^2`. Per la seconda, se `v<=4/3`,

\[
1+\frac v2-\left(\eta+\frac{v}{1+\eta^2}\right)
=(1-\eta)\left(1-\frac{v(1+\eta)}{2(1+\eta^2)}\right)\ge0,
\]

poiché

\[
2(1+\eta^2)-\frac43(1+\eta)
=\frac23(3\eta^2-2\eta+1)>0.
\]

Se `v>=4/3`, allora `(1+v/2)^2<=1+v^2`, quindi la prima
coordinata è sotto il secondo ramo di G. Per la seconda coordinata,
se `eta<=1/2`, Cauchy e

\[
\eta^2+(1+\eta^2)^{-2}\le1
\]

danno il bound. L'ultima disuguaglianza segue, con `z=eta^2<=1/4`,
da `z+(1+z)^(-2)-1=z(z^2+z-1)/(1+z)^2<=0`.
Se invece `eta>=1/2`, la derivata in `v` della differenza

\[
\sqrt{1+v^2}-\left(\eta+\frac{v}{1+\eta^2}\right)
\]

è almeno `4/5-1/(1+eta^2)>=0` per `v>=4/3`.
Al bordo `v=4/3`, tale differenza è non negativa per il primo caso,
e i due rami di G coincidono a `5/3`. Questo conclude.

Il bound logaritmico sul costo del minimo problematico è pertanto

\[
h_r(X)=
\begin{cases}
\frac12\ln(1+r^2/X^2),&0<X\le3r/4,\\
\ln(1+r/(2X)),&X\ge3r/4.
\end{cases} \tag{h}
\]

È continuo, positivo e decrescente. Per ogni `X>0` è strettamente
minore del precedente `ln(1+r/X)`. Il caso `a=0` corrisponde a
`eta=0` ed è compreso. Il cambio simultaneo di segno del vettore
non cambia alcun modulo.

## 3. Integrazione del termine lineare

La monotonia della restrizione e l'integrazione degli indicatori danno

\[
0\le\nu_{qp}(M)-\nu_q(M)
\le\int_0^\infty F_{\rm bad}(X)(-h_r'(X))\,dX.
\]

Il centro `-(1/2)ln M` si cancella fra i due marchi. L'eventuale
salto della derivata al raccordo non introduce un atomo, perché
`h_r` è continuo. I maggioranti sono integrabili ai due estremi.
Si separa `min(1,AX+Bsqrt(X))<=min(1,AX)+Bsqrt(X)`.

Ponendo `z=X/r`, il nucleo sul ramo inferiore è
`1/[z(1+z^2)]`; su quello superiore è `1/[z(2z+1)]`.
Per `y=rA<=4/3`, il punto `1/y` cade dopo `3/4`. Quindi

\[
\begin{split}
I_A&=\int_0^\infty\min(1,AX)(-h_r'(X))\,dX\\
&=y\arctan(3/4)
+\frac y2\ln\!\left(\frac{2+y}{(5/2)y}\right)
+\ln(1+y/2)\\
&=\Phi(y/2)+yK,
\qquad K=\arctan(3/4)-\frac12\ln(5/2).
\end{split}
\]

**CL3.** `verify_sharpened_pe.py` certifica esattamente

\[
0<K<0.185355807044945<\frac3{16}.
\]

Il bound di arctan usa la serie alternata fino all'indice pari 20:
la somma è un upper, e sottrarre il termine successivo dà un lower.
Il logaritmo è racchiuso con la serie razionale e il resto esplicito
di `GATES70.md`. Il decimale sopra è esterno e non decide il segno.

## 4. Integrazione del termine radicale

Si scrive `h0(X)=ln(1+r/(2X))` e `g=h_r-h0`.
Per `X>=3r/4` si ha `g=0`; sull'intervallo inferiore,

\[
g(X)=\ln\!\left(\frac{\sqrt{X^2+r^2}}{X+r/2}\right),
\qquad 0\le g(X)\le\ln2.
\]

Infatti `g` decresce dal limite `ln2` a zero: la sua derivata
è `r(X-2r)/[(X^2+r^2)(2X+r)]<0` in tale intervallo.
Il termine `h0` dà l'integrale esatto `pi*sqrt(r/2)`.
Per il termine restante, integrando per parti,

\[
\int_0^\infty\sqrt X(-g'(X))\,dX
=\frac12\int_0^{3r/4}\frac{g(X)}{\sqrt X}\,dX
\le\frac{\sqrt{3r}}2\ln2.
\]

I termini al bordo sono zero, incluso il limite `sqrt(X)*g(X)`
in zero. Pertanto

\[
\int_0^\infty\sqrt X(-h_r'(X))\,dX
\le\sqrt r\left(\frac\pi{\sqrt2}
+\frac{\sqrt3}{2}\ln2\right). \tag{J}
\]

### Controllo alternativo del resto

Una seconda decomposizione interna conferma il bound, con margine.
Sia `g0(X)=(1/2)ln(1+r^2/X^2)`. La differenza `D=h_r-g0`
è zero sotto `3r/4`, non negativa sopra, e soddisfa
`D(X)<=ln(1+r/(2X))<=r/(2X)`. Per parti,

\[
\int_0^\infty\sqrt X(-D'(X))\,dX
=\frac12\int_{3r/4}^{\infty}\frac{D(X)}{\sqrt X}\,dX
\le\frac{\sqrt r}{\sqrt3}.
\]

L'integrale del nucleo di `g0` è `pi*sqrt(r/2)`, tramite
`X=r t^2` e l'integrale elementare
`int_0^infinity dt/(1+t^4)=pi/(2sqrt2)`.
Quest'ultimo si ricava anche sommando la propria trasformata `t->1/t`
e usando le due quadratiche `t^2+sqrt2*t+1`,
`t^2-sqrt2*t+1`. Poiché `ln2>2/3`, segue
`1/sqrt3<(sqrt3/2)ln2`, e la seconda decomposizione sostiene J.
Non si usa questo margine per cambiare le costanti o i cutoff.

## 5. Lemma medio generale e coefficienti esatti

Dai §§3–4, sotto le ipotesi del §1,

\[
\nu_{qp}(M)-\nu_q(M)
\le\Phi(rA/2)+rA K
+B\sqrt r\left(\frac\pi{\sqrt2}
+\frac{\sqrt3}{2}\ln2\right).
\tag{PE-sharp}
\]

Per ottenere coefficienti razionali, il verificatore certifica
`sqrt2<99/70`, `sqrt3<26/15`, `ln2<25/36`; si conserva
`pi<22/7`, già dimostrato nell'integrale positivo della prova70.

Per `r=4`, il coefficiente radicale universale è

\[
\pi\sqrt2+\sqrt3\ln2
<\frac{22}{7}\frac{99}{70}+\frac{26}{15}\frac{25}{36}
=\frac{74731}{13230}.
\]

Per `r=6` è

\[
\pi\sqrt3+\frac{3\sqrt2}{2}\ln2
<\frac{22}{7}\frac{26}{15}
+\frac32\frac{99}{70}\frac{25}{36}
=\frac{11627}{1680}.
\]

Con P70 e `K<3/16` segue, per ogni `p>4` primo, `p` non
divisore di 70 e `p|M`,

\[
\boxed{\displaystyle
\nu_{70p}(M)-\nu_{70}(M)
<\Phi\!\left(\frac{12}{7p}\right)
+\frac9{14p}+\frac{971503}{44247\sqrt M}.}
\tag{PE70-sharp}
\]

Con P770, per ogni `ell>6` primo, `ell` non divisore di 770 e
`ell|M`, si ottiene

\[
\boxed{\displaystyle
\nu_{770\ell}(M)-\nu_{770}(M)
<\Phi\!\left(\frac{180}{77\ell}\right)
+\frac{135}{154\ell}+\frac{1755677}{43120\sqrt M}.}
\tag{PE770-sharp}
\]

In entrambi i casi `y<=4/3` segue già da `p>r`; non è richiesta
la divisibilità per 70 o 770 finché non si applica il relativo
ancoraggio marcato. La differenza di medie rimane non negativa
per monotonia della restrizione.

Per ogni nucleo `q=2uv`, con `u,v` primi dispari distinti,
la prova senza tabelle di `REPAIR70_PROOF.md`, §7, dà
`r=4`, `A=(5/2)rho/p`, `B=8/sqrt(M)`, `rho<1/2`.
Quindi, per `p>4`, `p` non divisore di `q` e `p|M`,

\[
\nu_{qp}(M)-\nu_q(M)
<\Phi\!\left(\frac5{2p}\right)
+\frac{15}{16p}+\frac{298924}{6615\sqrt M}.
\tag{PE2uv-sharp}
\]

Anche questo è un costo di inserimento: non fornisce da solo un
ancoraggio uniforme per `nu_q`.

## 6. Impiego sulle rotte già fissate

Sommare il medesimo ancoraggio marcato e dividere per **ln B**
produce un secondo bound per Delta. Con gli ancoraggi non negativi
`7/50`, `0`, `9/80`, il numeratore dell'espressione sufficiente è
positivo. Si possono usare minoranti razionali di `ln B` e
`sqrt(M)` senza ambiguità di segno. Se si sostituisse un ancoraggio
negativo più preciso, andrebbe gestito correttamente il segno del
numeratore nel rapporto di intervalli.

Il verificatore ricertifica soltanto i bordi già fissati:

| Rotta invariata | Nuovo upper decimale esterno per Delta |
|---|---:|
| 70p: p>=1000, tutti N | 0.021193904550578 |
| 70p: p>=73, N>=2 | 0.014633752502645 |
| 70p: p>=73, B>=3000000, tutti N | 0.018327399494273 |
| 770ell: ell>=300, tutti N | 0.019910491075429 |

Gli upper e i confronti con le precedenti espressioni sono frazioni
esatte nel JSON; i decimali non sono valori effettivi di Delta.
La monotonia in `p,M,B` dell'espressione nuova segue dalla positività
e monotonia di Phi e degli altri termini. La stessa formula può
essere valutata sulle medesime basi del complemento finito, senza
cambiare il suo inventario o rigenerare profili geometrici.

`verify_sharpened_pe.py` espone le funzioni numeriche
`sharp_pe70_upper(p,M)` e `sharp_pe770_upper(p,M)`, che restituiscono
un upper razionale del solo costo di inserimento, con radice minorata
da `isqrt(M)`. La divisione per `ln B` viene fatta dal chiamante.
Le funzioni accettano anche parametri scalari di bordo; l'enunciato
matematico sulle medie conserva le ipotesi di primalità e divisibilità.

```bash
python verify_sharpened_pe.py
```

I 48 controlli puntuali a scale miste sono regressioni finite del
bound geometrico. La validità per ogni `eta,v` è dimostrata nel §2,
non inferita da quei controlli. I gate accettati in precedenza,
gli output dei profili e gli inventari rimangono invariati.
